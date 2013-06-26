#ifndef BAM_MATE_ITER_H
#define BAM_MATE_ITER_H

#include "samtools/sam.h"

typedef struct _bam_mate_iter_t *bam_mate_iter_t;

enum MATE_RESOLUTION { SINGLE_END, MATES_FOUND, MATE_UNAVAILABLE, MATES_NOT_FOUND };

typedef struct {
    int n;
    bam1_t **bams;
} bam_mates_t;

#ifdef __cplusplus
extern "C" {
#endif

bam_mates_t *bam_mates_new();
void bam_mates_realloc(bam_mates_t *mates, int n);
void bam_mates_destroy(bam_mates_t *mates);

int samread_mate(samfile_t *fp, bam_mates_t *b);
int bam_mate_fetch(bamFile fb, const bam_index_t *idx, int tid, int beg, 
                   int end, void *data, bam_fetch_f func);
bam_mate_iter_t bam_mate_query(const bam_index_t *, int tid, int beg, int end);
int bam_mate_read(bamFile fp, bam_mate_iter_t iter, bam_mates_t *mates);
void bam_mate_destroy(bam_mate_iter_t iter);

#ifdef __cplusplus
}
#endif

#endif
