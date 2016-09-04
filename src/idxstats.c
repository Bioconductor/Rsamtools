#include "bamfile.h"
#include "utilities.h"          /* _checkext */
#include "samtools/khash.h"
#include "samtools/ksort.h"

extern SEXP BAMFILE_TAG;

/* 
   The struct and macro definitions are from samtools/bamindex.c and
   are a hack exposing internal business.

 */

#define BAM_MAX_BIN 37450 // =(8^6-1)/7+1

typedef struct {
	uint64_t u, v;
} pair64_t;

#define pair64_lt(a,b) ((a).u < (b).u)

typedef struct {
	uint32_t m, n;
	pair64_t *list;
} bam_binlist_t;

typedef struct {
	int32_t n, m;
	uint64_t *offset;
} bam_lidx_t;

KHASH_MAP_INIT_INT(idxstat_hash, bam_binlist_t)

struct __bam_index_t {
	int32_t n;
	uint64_t n_no_coor; // unmapped reads without coordinate
	khash_t(idxstat_hash) **index;
	bam_lidx_t *index2;
};

SEXP idxstats_bamfile(SEXP ext)
{
    BAM_FILE bfile;
    bamFile fp;
    bam_header_t *header;
    bam_index_t *idx;
    SEXP result, name, len, map, unmap;

    _checkext(ext, BAMFILE_TAG, "idxstats");
    bfile = BAMFILE(ext);
    fp = bfile->file->x.bam;
    bam_seek(fp, 0, 0);
    header = bam_header_read(fp);
    idx = bfile->index;

    result = PROTECT(Rf_allocVector(VECSXP, 4));
    name = Rf_allocVector(STRSXP, idx->n); SET_VECTOR_ELT(result, 0, name);
    len  = Rf_allocVector(INTSXP, idx->n); SET_VECTOR_ELT(result, 1, len);
    map = Rf_allocVector(REALSXP, idx->n); SET_VECTOR_ELT(result, 2, map);
    unmap = Rf_allocVector(REALSXP, idx->n); SET_VECTOR_ELT(result, 3, unmap);
    
    for (int i = 0; i < idx->n; ++i) {
        khint_t k;
        khash_t(idxstat_hash) *h = idx->index[i];
        SET_STRING_ELT(name, i, mkChar(header->target_name[i]));
        INTEGER(len)[i] = header->target_len[i];
        k = kh_get(idxstat_hash, h, BAM_MAX_BIN);
        if (k != kh_end(h)) {
            REAL(map)[i] = (double) kh_val(h, k).list[1].u;
            REAL(unmap)[i] = (double) kh_val(h, k).list[1].v;
        } else {
            REAL(map)[i] = REAL(unmap)[i] = 0;
        }
    }

    UNPROTECT(1);
    return result;
}
