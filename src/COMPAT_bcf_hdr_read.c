#include "COMPAT_bcf_hdr_read.h"

#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/hfile.h>
#include <htslib/bgzf.h>
#include <errno.h>

/* WARNING: Between old samtools (version < 0.1.19) and htslib 1.7, the
   behavior of vcf_hdr_read() has changed as follow: the new version now
   checks for the presence of an existing tabix index and uses it to add
   contigs not listed in the header.
   COMPAT_vcf_hdr_read() below is a re-implementation of the old vcf_hdr_read()
   i.e. it does NOT try to use a possible existing tabix index in order to
   add missing contigs.
   We use it in by scan_bcf_header() (instead of vcf_hdr_read()) so its
   behavior remains the same as before the migration of Rsamtools to
   htslib 1.7 (in particular the Reference component of the returned list
   remains the same).

   Testing:
     library(Rsamtools)
     bcf <- system.file("extdata", "ex1.bcf.gz", package="Rsamtools")
     scanBcfHeader(bcf)

     library(VariantAnnotation)
     vcf <- system.file("extdata", "structural.vcf",
                        package="VariantAnnotation")
     bgz <- bgzip(vcf, tempfile())
     target <- scanBcfHeader(BcfFile(vcf, character(0)))[[1]]
     stopifnot(identical(target,
                         scanBcfHeader(BcfFile(bgz, character(0)))[[1]]))
     idx <- indexTabix(bgz, "vcf")
     stopifnot(identical(target,
                         scanBcfHeader(BcfFile(bgz, character(0)))[[1]]))
     stopifnot(identical(target,
                         scanBcfHeader(BcfFile(bgz, idx))[[1]]))
*/
static bcf_hdr_t *COMPAT_vcf_hdr_read(htsFile *fp)
{
    kstring_t txt, *s = &fp->line;
    int ret;
    bcf_hdr_t *h;
    h = bcf_hdr_init("r");
    if (!h) {
        hts_log_error("Failed to allocate bcf header");
        return NULL;
    }
    txt.l = txt.m = 0; txt.s = 0;
    while ((ret = hts_getline(fp, KS_SEP_LINE, s)) >= 0) {
        if (s->l == 0) continue;
        if (s->s[0] != '#') {
            hts_log_error("No sample line");
            goto error;
        }
        if (s->s[1] != '#' && fp->fn_aux) { // insert contigs here
            kstring_t tmp = { 0, 0, NULL };
            hFILE *f = hopen(fp->fn_aux, "r");
            if (f == NULL) {
                hts_log_error("Couldn't open \"%s\"", fp->fn_aux);
                goto error;
            }
            while (tmp.l = 0, kgetline(&tmp, (kgets_func *) hgets, f) >= 0) {
                char *tab = strchr(tmp.s, '\t');
                if (tab == NULL) continue;
                kputs("##contig=<ID=", &txt); kputsn(tmp.s, tab - tmp.s, &txt);
                kputs(",length=", &txt); kputl(atol(tab), &txt);
                kputsn(">\n", 2, &txt);
            }
            free(tmp.s);
            if (hclose(f) != 0) {
                hts_log_warning("Failed to close %s", fp->fn_aux);
            }
        }
        kputsn(s->s, s->l, &txt);
        kputc('\n', &txt);
        if (s->s[1] != '#') break;
    }
    if ( ret < -1 ) goto error;
    if ( !txt.s )
    {
        hts_log_error("Could not read the header");
        goto error;
    }
    if ( bcf_hdr_parse(h, txt.s) < 0 ) goto error;
    free(txt.s);
    return h;

 error:
    free(txt.s);
    if (h) bcf_hdr_destroy(h);
    return NULL;
}

bcf_hdr_t *COMPAT_bcf_hdr_read(htsFile *hfp)
{
    if (hfp->format.format == vcf)
        return COMPAT_vcf_hdr_read(hfp);
    if (hfp->format.format != bcf) {
        hts_log_error("Input is not detected as bcf or vcf format");
        return NULL;
    }

    assert(hfp->is_bgzf);

    BGZF *fp = hfp->fp.bgzf;
    uint8_t magic[5];
    bcf_hdr_t *h;
    h = bcf_hdr_init("r");
    if (!h) {
        hts_log_error("Failed to allocate bcf header");
        return NULL;
    }
    if (bgzf_read(fp, magic, 5) != 5)
    {
        hts_log_error("Failed to read the header (reading BCF in text mode?)");
        bcf_hdr_destroy(h);
        return NULL;
    }
    if (strncmp((char*)magic, "BCF\2\2", 5) != 0)
    {
        if (!strncmp((char*)magic, "BCF", 3))
            hts_log_error("Invalid BCF2 magic string: only BCFv2.2 is supported");
        else
            hts_log_error("Invalid BCF2 magic string");
        bcf_hdr_destroy(h);
        return NULL;
    }
    uint8_t buf[4];
    size_t hlen;
    char *htxt = NULL;
    if (bgzf_read(fp, buf, 4) != 4) goto fail;
    hlen = buf[0] | (buf[1] << 8) | (buf[2] << 16) | (buf[3] << 24);
    if (hlen >= SIZE_MAX) { errno = ENOMEM; goto fail; }
    htxt = (char*)malloc(hlen + 1);
    if (!htxt) goto fail;
    if (bgzf_read(fp, htxt, hlen) != hlen) goto fail;
    htxt[hlen] = '\0'; // Ensure htxt is terminated
    if ( bcf_hdr_parse(h, htxt) < 0 ) goto fail;
    free(htxt);
    return h;
 fail:
    hts_log_error("Failed to read BCF header");
    free(htxt);
    bcf_hdr_destroy(h);
    return NULL;
}

