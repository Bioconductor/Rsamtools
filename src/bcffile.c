#include <string.h>
#include <errno.h>
#include "samtools/kstring.h"
#include "bcftools/bcf.h"
#include "bcffile.h"
#include "utilities.h"

#define smkChar(x) ((x) ? mkChar(x) : R_NaString)

struct typemap {
    uint32_t *id;
    int *type, *mult;
    int len;
};

enum {
    BCF_HDR_REF = 0, BCF_HDR_SAMPLE, BCF_HDR_HEADER, BCF_HDR_LAST
};

static const char *BCF_HDR_NM[] = { "Reference", "Sample", "Header" };

enum {
    BCF_TID = 0, BCF_POS, BCF_ID, BCF_REF, BCF_ALT, BCF_QUAL,
    BCF_FLT, BCF_INFO, BCF_FMT, BCF_GENO, BCF_RECS_PER_RANGE,
    BCF_LAST
};

enum {
    BCF_TYPE_Integer = 1, BCF_TYPE_Float, BCF_TYPE_Character,
    BCF_TYPE_String, BCF_TYPE_Flag, BCF_TYPE_Last
};

static SEXP BCFFILE_TAG = NULL;
static const int BCF_BUFSIZE_GROW = 100000;	/* initial # records */

static bcf_t *_bcf_tryopen(const char *fname, const char *mode)
{
    return vcf_open(fname, mode);
}

static bcf_idx_t *_bcf_idx_load(const char *fname)
{
    return bcf_idx_load(fname);
}

static void _bcf_close(bcf_t * bcf, Rboolean errmsg)
{
    int err = vcf_close(bcf);
    if (0 != err && errmsg) {
        if (Z_ERRNO == err) {
            err = errno;
            Rf_error("_bcf_close file system error (%d): %s",
                     err, strerror(err));
        }
        Rf_error("_bcf_close error (%d)", err);
    }
}

static void _bcffile_close(SEXP ext)
{
    _BCF_FILE *bfile = BCFFILE(ext);
    if (NULL != bfile->file)
        vcf_close(bfile->file);
    if (NULL != bfile->index)
        bcf_idx_destroy(bfile->index);
    bfile->file = NULL;
    bfile->index = NULL;
}

static void _bcffile_finalizer(SEXP ext)
{
    if (NULL == R_ExternalPtrAddr(ext))
        return;
    _bcffile_close(ext);
    _BCF_FILE *bfile = BCFFILE(ext);
    Free(bfile);
    R_SetExternalPtrAddr(ext, NULL);
}

SEXP bcffile_init()
{
    BCFFILE_TAG = install("BcfFile");
    return R_NilValue;
}

SEXP bcffile_open(SEXP filename, SEXP indexname, SEXP filemode)
{
    _scan_checknames(filename, indexname, filemode);

    _BCF_FILE *bfile = Calloc(1, _BCF_FILE);

    bfile->file = NULL;
    if (0 != Rf_length(filename)) {
        const char *cfile = translateChar(STRING_ELT(filename, 0));
        bfile->file = _bcf_tryopen(cfile, CHAR(STRING_ELT(filemode, 0)));
        if (NULL == bfile->file) {
            Free(bfile);
            Rf_error("'open' BCF failed\n  filename: %s", cfile);
        }
    }

    bfile->index = NULL;
    if (0 != Rf_length(indexname) && !bfile->file->is_vcf) {
        const char *cindex = translateChar(STRING_ELT(indexname, 0));
        bfile->index = _bcf_idx_load(cindex);
        if (NULL == bfile->index) {
            _bcf_close(bfile->file, FALSE);
            Free(bfile);
            Rf_error("'open' BCF index failed\n  indexname: %s\n", cindex);
        }
    }

    SEXP ext = PROTECT(R_MakeExternalPtr(bfile, BCFFILE_TAG, filename));
    R_RegisterCFinalizerEx(ext, _bcffile_finalizer, TRUE);
    UNPROTECT(1);

    return ext;
}

SEXP bcffile_close(SEXP ext)
{
    _scan_checkext(ext, BCFFILE_TAG, "close");
    _bcffile_close(ext);
    return ext;
}

SEXP bcffile_isopen(SEXP ext)
{
    SEXP ans = ScalarLogical(FALSE);
    if (NULL != BCFFILE(ext)) {
        _scan_checkext(ext, BCFFILE_TAG, "isOpen");
        if (BCFFILE(ext)->file)
            ans = ScalarLogical(TRUE);
    }
    return ans;
}

SEXP bcffile_isvcf(SEXP ext)
{
    SEXP ans = ScalarLogical(FALSE);
    if (NULL != BCFFILE(ext)) {
        _scan_checkext(ext, BCFFILE_TAG, "isVcf");
        if (BCFFILE(ext)->file && BCFFILE(ext)->file->is_vcf)
            ans = ScalarLogical(TRUE);
    }
    return ans;
}

/* implementation */
static int _bcf_ans_grow(SEXP ans, R_len_t sz, int n_smpl)
{
    R_len_t n = sz;
    if (0 <= sz)
        n += Rf_length(VECTOR_ELT(ans, BCF_TID));
    else
        n *= -1;
    for (int i = 0; i < BCF_LAST; ++i) {
        SEXP elt = VECTOR_ELT(ans, i);
        switch (i) {
        case BCF_GENO:
            for (int i = 0; i < Rf_length(elt); ++i) {
                SEXP g = VECTOR_ELT(elt, i);
                SEXP dim = GET_DIM(g);
                if (R_NilValue == dim) {
                    g = Rf_lengthgets(g, n);
                    SET_VECTOR_ELT(elt, i, g);	/* protect */
                } else {
                    PROTECT(dim);
                    g = Rf_lengthgets(g, n_smpl * n);
                    SET_VECTOR_ELT(elt, i, g);	/* protect */
                    INTEGER(dim)[0] = n_smpl;
                    INTEGER(dim)[1] = n;
                    SET_DIM(g, dim);
                    UNPROTECT(1);
                }
            }
            break;
        case BCF_RECS_PER_RANGE:
            break;
        default:
            elt = Rf_lengthgets(elt, n);
            SET_VECTOR_ELT(ans, i, elt);	/* protect */
            break;
        }
    }
    return n;
}

static int _bcf_sync1(bcf1_t * b)
{
    /* called when no FORMAT / GENO fields present;
       from bcftools/bcf.c:bcf_sync */
    char *p, *tmp[5];
    int n;
    for (p = b->str, n = 0; p < b->str + b->l_str; ++p) {
        if (*p == 0 && p + 1 != b->str + b->l_str) {
            if (n == 5) {
                ++n;
                break;
            } else
                tmp[n++] = p + 1;
        }
    }
    if (n != 4)
        return -1;
    b->ref = tmp[0];
    b->alt = tmp[1];
    b->flt = tmp[2];
    b->info = tmp[3];
    b->fmt = 0;
    if (*b->alt == 0)
        b->n_alleles = 1;
    else {
        for (p = b->alt, n = 1; *p; ++p)
            if (*p == ',')
                ++n;
        b->n_alleles = n + 1;
    }
    b->n_gi = 0;
    return 0;
}

static void _bcf_gi2sxp(SEXP geno, const int i_rec, const bcf_hdr_t * h,
                        bcf1_t * b)
{
    SEXP nm = GET_NAMES(geno);
    /* FIXME: more flexible geno not supported by bcftools */

    if (b->n_gi == 0)
        return;

    /* from bcftools/bcf.c */
    for (int i = 0; i < b->n_gi; ++i) {
        const int off = i_rec * h->n_smpl;
        SEXP g;
        int t;

        for (t = 0; t < Rf_length(nm); ++t)
            if (bcf_str2int(CHAR(STRING_ELT(nm, t)), 2) == b->gi[i].fmt)
                break;
        if (Rf_length(nm) <= t)
            Rf_error("failed to find fmt encoded as '%d'", b->gi[i].fmt);
        g = VECTOR_ELT(geno, t);

        if (b->gi[i].fmt == bcf_str2int("PL", 2)) {
            const int x = b->n_alleles * (b->n_alleles + 1) / 2;
            SEXP pl = Rf_allocMatrix(INTSXP, x, h->n_smpl);
            SET_VECTOR_ELT(g, off, pl);	/* protect */
            for (int j = 0; j < h->n_smpl; ++j) {
                uint8_t *d = (uint8_t *) b->gi[i].data + j * x;
                for (int k = 0; k < x; ++k)
                    INTEGER(pl)[j * x + k] = d[k];
            }
        } else if (b->gi[i].fmt == bcf_str2int("DP", 2)) {
            int *dp = INTEGER(g) + off;
            for (int j = 0; j < h->n_smpl; ++j)
                *dp++ = ((uint16_t *) b->gi[i].data)[j];
        } else if (b->gi[i].fmt == bcf_str2int("GQ", 2) ||
                   b->gi[i].fmt == bcf_str2int("SP", 2)) {
            int *gq = INTEGER(g) + off;
            for (int j = 0; j < h->n_smpl; ++j)
                *gq++ = ((uint8_t *) b->gi[i].data)[j];
        } else if (b->gi[i].fmt == bcf_str2int("GT", 2)) {
            int idx = off;
            char s[4];
            s[3] = '\0';
            for (int j = 0; j < h->n_smpl; ++j) {
                int y = ((uint8_t *) b->gi[i].data)[j];
                if (y >> 7 & 1)
                    SET_STRING_ELT(g, idx++, mkChar("./."));
                else {
                    s[0] = '0' + (y >> 3 & 7);
                    s[1] = "/|"[y >> 6 & 1];
                    s[2] = '0' + (y & 7);
                    SET_STRING_ELT(g, idx++, mkChar(s));
                }
            }
        } else if (b->gi[i].fmt == bcf_str2int("GL", 2)) {
            const int x = b->n_alleles * (b->n_alleles + 1) / 2;
            SEXP gl = Rf_allocMatrix(REALSXP, x, h->n_smpl);
            SET_VECTOR_ELT(g, off, gl);	/* protect */
            for (int j = 0; j < h->n_smpl; ++j) {
                float *d = (float *) b->gi[i].data + j * x;
                for (int k = 0; k < x; ++k)
                    REAL(gl)[j * x + k] = d[k];
            }
        }
    }
}

SEXP scan_bcf_header(SEXP ext)
{
    _scan_checkext(ext, BCFFILE_TAG, "scanBcfHeader");
    bcf_t *bcf = BCFFILE(ext)->file;
    if (!bcf->is_vcf && 0 != bgzf_seek(bcf->fp, 0, SEEK_SET))
        Rf_error("internal: failed to 'seek'");
    bcf_hdr_t *hdr = vcf_hdr_read(bcf);

    SEXP ans = PROTECT(NEW_LIST(BCF_HDR_LAST));
    SET_VECTOR_ELT(ans, BCF_HDR_REF, NEW_STRING(hdr->n_ref));
    SET_VECTOR_ELT(ans, BCF_HDR_SAMPLE, NEW_STRING(hdr->n_smpl));
    /* count header text lines */
    const char *c;
    char *s;
    int n_hdr = 0;
    if (hdr->l_txt)
        for (c = hdr->txt; *c != '\0'; ++c)
            if (*c == '\n')
                ++n_hdr;
    SET_VECTOR_ELT(ans, BCF_HDR_HEADER, NEW_STRING(n_hdr));

    int i;
    SEXP x = VECTOR_ELT(ans, BCF_HDR_REF);
    for (i = 0; i < hdr->n_ref; ++i)
        SET_STRING_ELT(x, i, mkChar(_rtrim(hdr->ns[i])));
    x = VECTOR_ELT(ans, BCF_HDR_SAMPLE);
    for (i = 0; i < hdr->n_smpl; ++i)
        SET_STRING_ELT(x, i, mkChar(_rtrim(hdr->sns[i])));
    x = VECTOR_ELT(ans, BCF_HDR_HEADER);
    if (hdr->l_txt) {
        char *txt = (char *) R_alloc(hdr->l_txt, sizeof(char));
        strncpy(txt, hdr->txt, hdr->l_txt);
        s = strtok(txt, "\n");
        for (i = 0; i < n_hdr; ++i) {
            SET_STRING_ELT(x, i, mkChar(_rtrim(s)));
            s = strtok(NULL, "\n");
        }
    }

    SEXP nm = NEW_CHARACTER(3);
    SET_NAMES(ans, nm);         /* protect */
    for (i = 0; i < BCF_HDR_LAST; ++i)
        SET_STRING_ELT(nm, i, mkChar(BCF_HDR_NM[i]));

    bcf_hdr_destroy(hdr);
    UNPROTECT(1);
    return ans;
}

int scan_bcf_range(bcf_t * bcf, bcf_hdr_t * hdr, SEXP ans, int tid, int start,
                   int end, int n)
{
    const int TID_BUFSZ = 8;
    static char *buf = NULL;
    bcf1_t *bcf1 = calloc(1, sizeof(bcf1_t));	/* free'd in bcf_destroy */
    if (NULL == bcf1)
        Rf_error("scan_bcf_region: failed to allocate memory");
    int sz = Rf_length(VECTOR_ELT(ans, BCF_TID));
    int res;
    if (NULL == buf)
        buf = Calloc(TID_BUFSZ, char);	/* leaks, but oh well */
    while (0 <= (res = vcf_read(bcf, hdr, bcf1))) {
        if (tid >= 0) {
            int pos = strlen(bcf1->ref);
            pos = bcf1->pos + (pos > 0 ? pos : 1);
            if (bcf1->tid != tid || bcf1->pos > end)
                break;
            if (!(pos >= start && end > bcf1->pos))
                continue;
        }
        if (n >= sz)
            sz = _bcf_ans_grow(ans, BCF_BUFSIZE_GROW, hdr->n_smpl);
        if (n >= sz) {
            bcf_destroy(bcf1);
            Rf_error("bcf_scan: failed to increase size; out of memory?");
        }
        if (hdr->ns)
            SET_STRING_ELT(VECTOR_ELT(ans, BCF_TID), n,
                           smkChar(hdr->ns[bcf1->tid]));
        else {
            snprintf(buf, TID_BUFSZ, "%d", bcf1->tid);
            SET_STRING_ELT(VECTOR_ELT(ans, BCF_TID), n, smkChar(buf));
        }
        if (bcf->is_vcf && NULL == bcf1->ref)
            if (_bcf_sync1(bcf1)) {
                bcf_destroy(bcf1);
                Rf_error("bcf_scan: unexpected number of fields in line %d",
                         n + 1);
            }
        INTEGER(VECTOR_ELT(ans, BCF_POS))[n] = bcf1->pos + 1;
        REAL(VECTOR_ELT(ans, BCF_QUAL))[n] = bcf1->qual;
        SET_STRING_ELT(VECTOR_ELT(ans, BCF_ID), n, smkChar(bcf1->str));
        SET_STRING_ELT(VECTOR_ELT(ans, BCF_REF), n, smkChar(bcf1->ref));
        SET_STRING_ELT(VECTOR_ELT(ans, BCF_ALT), n, smkChar(bcf1->alt));
        SET_STRING_ELT(VECTOR_ELT(ans, BCF_FLT), n, smkChar(bcf1->flt));
        SET_STRING_ELT(VECTOR_ELT(ans, BCF_INFO), n, smkChar(bcf1->info));
        SET_STRING_ELT(VECTOR_ELT(ans, BCF_FMT), n, smkChar(bcf1->fmt));
        _bcf_gi2sxp(VECTOR_ELT(ans, BCF_GENO), n, hdr, bcf1);
        if (bcf->is_vcf)
            bcf1->ref = NULL;
        ++n;
    }
    bcf_destroy(bcf1);
    return n;
}

SEXP scan_bcf(SEXP ext, SEXP space, SEXP tmpl)
{
    _scan_checkparams(space, R_NilValue, R_NilValue);
    _scan_checkext(ext, BCFFILE_TAG, "scanBcf");
    bcf_t *bcf = BCFFILE(ext)->file;
    bcf_idx_t *idx = BCFFILE(ext)->index;
    if (!bcf->is_vcf && 0 != bgzf_seek(bcf->fp, 0, SEEK_SET))
        Rf_error("internal: failed to 'seek' on bcf file");
    bcf_hdr_t *hdr = vcf_hdr_read(bcf);
    if (NULL == hdr)
        Rf_error("failed to read header; wrong 'mode' or corrupt file?");

    int n = 0;
    tmpl = PROTECT(Rf_duplicate(tmpl));

    if (R_NilValue == space) {
        SET_VECTOR_ELT(tmpl, BCF_RECS_PER_RANGE, NEW_INTEGER(1));
        n = scan_bcf_range(bcf, hdr, tmpl, -1, -1, -1, n);
        INTEGER(VECTOR_ELT(tmpl, BCF_RECS_PER_RANGE))[0] = n;
    } else {
        SEXP spc = VECTOR_ELT(space, 0);
        const int
        *start = INTEGER(VECTOR_ELT(space, 1)),
            *end = INTEGER(VECTOR_ELT(space, 2)), nspc = Rf_length(spc);
        void *str2id = bcf_build_refhash(hdr);
        SEXP nrec = NEW_INTEGER(nspc);
        SET_VECTOR_ELT(tmpl, BCF_RECS_PER_RANGE, nrec);

        for (int i = 0; i < nspc; ++i) {
            int tid = bcf_str2id(str2id, CHAR(STRING_ELT(spc, i)));
            if (tid < 0) {
                bcf_str2id_destroy(str2id);
                Rf_error("'space' not in file: %s", CHAR(STRING_ELT(spc, i)));
            }
            uint64_t off = bcf_idx_query(idx, tid, start[i]);
            if (off == 0) {
                INTEGER(nrec)[i] = 0;
                continue;
            }
            bgzf_seek(bcf->fp, off, SEEK_SET);
            n = scan_bcf_range(bcf, hdr, tmpl, tid, start[i], end[i], n);
            if (i == 0)
                INTEGER(nrec)[i] = n;
            else
                INTEGER(nrec)[i] = n - INTEGER(nrec)[i - 1];
        }
        bcf_str2id_destroy(str2id);
    }
    _bcf_ans_grow(tmpl, -1 * n, hdr->n_smpl);

    UNPROTECT(1);
    return tmpl;
}

int _as_bcf(bcf_t * fin, const char *dict, bcf_t * fout)
{
    bcf1_t *b = calloc(1, sizeof(bcf1_t));	/* free'd in bcf_destroy */
    if (NULL == b)
        Rf_error("_as_bcf: failed to allocate memory");
    bcf_hdr_t *hin, *hout;
    int r, count = 0;

    hin = hout = vcf_hdr_read(fin);
    vcf_dictread(fin, hin, dict);
    vcf_hdr_write(fout, hout);
    while (0 <= (r = vcf_read(fin, hin, b))) {
        if (NULL == b->ref)
            Rf_error("cannot (yet) coerce VCF files without FORMAT");
        vcf_write(fout, hout, b);
        count++;
    }

    if (hin != hout)
        bcf_hdr_destroy(hout);
    bcf_hdr_destroy(hin);
    bcf_destroy(b);

    return r >= -1 ? count : -1 * count;
}

SEXP as_bcf(SEXP file, SEXP dictionary, SEXP destination)
{
    if (!IS_CHARACTER(file) || 1 != LENGTH(file))
        Rf_error("'file' must be character(1)");
    if (!IS_CHARACTER(dictionary) || 1 != LENGTH(dictionary))
        Rf_error("'dictionary' must be character(1)");
    if (!IS_CHARACTER(destination) || 1 != LENGTH(destination))
        Rf_error("'destination' must be character(1)");

    bcf_t *fin = _bcf_tryopen(translateChar(STRING_ELT(file, 0)), "r");
    if (NULL == fin)
        Rf_error("failed to open VCF 'file'");

    bcf_t *fout = _bcf_tryopen(translateChar(STRING_ELT(destination, 0)), "wb");
    if (NULL == fout)
        Rf_error("failed to open BCF 'destination'");

    int count = _as_bcf(fin, translateChar(STRING_ELT(dictionary, 0)), fout);

    _bcf_close(fin, FALSE);
    _bcf_close(fout, FALSE);
    if (count < 0)
        Rf_error("truncated input file at record %d", -1 * count);

    return destination;
}

SEXP index_bcf(SEXP file)
{
    if (!IS_CHARACTER(file) || 1 != LENGTH(file))
        Rf_error("'file' must be character(1)");
    const char *fbcf = translateChar(STRING_ELT(file, 0));
    int status = bcf_idx_build(fbcf);
    if (0 != status)
        Rf_error("failed to build index");
    char *fidx = (char *) R_alloc(strlen(fbcf) + 5, sizeof(char));
    sprintf(fidx, "%s.bci", fbcf);
    return mkString(fidx);
}
