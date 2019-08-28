#include "bamfile.h"
#include "utilities.h"          /* _checkext */
#include <htslib/khash.h>
#include <htslib/ksort.h>
#include "hts_utilities.h"

extern SEXP BAMFILE_TAG;

/*
 * Adapted from samtools-1.9/bam_index.c:106 - 165
 *
 * Cram indices do not contain mapped/unmapped record counts, so we have to
 * decode each record and count.  However we can speed this up as much as
 * possible by using the required fields parameter.
 */

int _idxstats_slow(htsFile *fp, bam_hdr_t *header, uint64_t (*counts)[2])
{
    int ret, last_tid = -2;
    bam1_t *b = bam_init1();

    if (hts_set_opt(fp, CRAM_OPT_REQUIRED_FIELDS, SAM_RNAME | SAM_FLAG))
        return -1;

    while ((ret = sam_read1(fp, header, b)) >= 0) {
        if (b->core.tid >= header->n_targets || b->core.tid < -1) {
            bam_destroy1(b);
            return -1;
        }

        if (b->core.tid != last_tid) {
            if (last_tid >= -1) {
                if (counts[b->core.tid][0] + counts[b->core.tid][1]) {
                    bam_destroy1(b);
                    hts_log_error("file is not position sorted");
                    return -1;
                }
            }
            last_tid = b->core.tid;
        }

        counts[b->core.tid][(b->core.flag & BAM_FUNMAP) ? 1 : 0]++;
    }

    bam_destroy1(b);

    if (hts_set_opt(fp, CRAM_OPT_REQUIRED_FIELDS, INT_MAX))
        return -1;

    return (ret == -1) ? 0 : -1;
}

int _idxstats_fast(bam_hdr_t *header, bam_index_t *idx, uint64_t (*counts)[2])
{
    uint32_t n = header->n_targets;

    for (int i = 0; i < n; ++i) {
        uint64_t mapped, unmapped;
        hts_idx_get_stat(idx, i, &mapped, &unmapped);
        counts[i][0] += mapped;
        counts[i][1] += unmapped;
    }
    counts[-1][0] = counts[-1][1] = 0;

    return 0;
}

SEXP idxstats_bamfile(SEXP ext)
{
    _checkext(ext, BAMFILE_TAG, "idxstats");
    BAM_FILE bfile = BAMFILE(ext);
    bam_hdr_t *header = bfile->header;

    int32_t n = header->n_targets;
    uint64_t (*count0)[2] = calloc(n + 1, sizeof(*count0));
    uint64_t (*counts)[2] = count0 + 1;

    int status;
    if (hts_get_format(bfile->file)->format != bam)
        status = _idxstats_slow(bfile->file, header, counts);
    else
        status = _idxstats_fast(header, bfile->index, counts);

    if (status < 0) {
        free(count0);
        Rf_error("'idxstats' failed");
    }

    SEXP result, name, len, map, unmap;

    result = PROTECT(Rf_allocVector(VECSXP, 4));
    name = Rf_allocVector(STRSXP, n + 1L); SET_VECTOR_ELT(result, 0, name);
    len  = Rf_allocVector(INTSXP, n + 1L); SET_VECTOR_ELT(result, 1, len);
    map = Rf_allocVector(REALSXP, n + 1L); SET_VECTOR_ELT(result, 2, map);
    unmap = Rf_allocVector(REALSXP, n + 1L); SET_VECTOR_ELT(result, 3, unmap);
    
    for (int i = 0; i < n; ++i) {
        SET_STRING_ELT(name, i, mkChar(header->target_name[i]));
        INTEGER(len)[i] = header->target_len[i];
        REAL(map)[i] = counts[i][0];
        REAL(unmap)[i] = counts[i][1];
    }
    /* unmapped reads */
    SET_STRING_ELT(name, n , mkChar("*"));
    INTEGER(len)[n] = 0;
    REAL(map)[n] = 0;
    REAL(unmap)[n] = hts_idx_get_n_no_coor(idx);

    free(count0);
    UNPROTECT(1);

    return result;
}
