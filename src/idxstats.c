#include "bamfile.h"
#include "utilities.h"          /* _checkext */
#include <htslib/khash.h>
#include <htslib/ksort.h>

extern SEXP BAMFILE_TAG;

SEXP idxstats_bamfile(SEXP ext)
{
    BAM_FILE bfile;
    bamFile fp;
    bam_header_t *header;
    bam_index_t *idx;
    int32_t n;
    SEXP result, name, len, map, unmap;

    _checkext(ext, BAMFILE_TAG, "idxstats");
    bfile = BAMFILE(ext);
    fp = bfile->file->x.bam;
    bam_seek(fp, 0, 0);
    header = bam_header_read(fp);
    idx = bfile->index;
    n = hts_idx_get_n(idx);

    result = PROTECT(Rf_allocVector(VECSXP, 4));
    name = Rf_allocVector(STRSXP, n); SET_VECTOR_ELT(result, 0, name);
    len  = Rf_allocVector(INTSXP, n); SET_VECTOR_ELT(result, 1, len);
    map = Rf_allocVector(REALSXP, n); SET_VECTOR_ELT(result, 2, map);
    unmap = Rf_allocVector(REALSXP, n); SET_VECTOR_ELT(result, 3, unmap);
    
    for (int i = 0; i < n; ++i) {
        uint64_t mapped, unmapped;
        SET_STRING_ELT(name, i, mkChar(header->target_name[i]));
        INTEGER(len)[i] = header->target_len[i];
        hts_idx_get_stat(idx, i, &mapped, &unmapped);
        REAL(map)[i] = (double) mapped;
        REAL(unmap)[i] = (double) unmapped;
    }

    UNPROTECT(1);
    return result;
}
