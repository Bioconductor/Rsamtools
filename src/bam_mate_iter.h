#ifndef BAM_MATE_ITER_H
#define BAM_MATE_ITER_H

#include "samtools/sam.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _bam_mate_iter_t *bam_mate_iter_t;

enum MATE_RESOLUTION { SINGLE_END, MATES_FOUND, MATE_UNAVAILABLE, MATES_NOT_FOUND };

typedef struct {
    int n;
    int mates;
    bam1_t **bams;
} bam_mates_t;

typedef void (*set_partition_f)(int partition, void *data);
typedef void (*set_mates_f)(int mate, void *data);

bam_mates_t *bam_mates_new();
void bam_mates_realloc(bam_mates_t *mates, int n);
void bam_mates_destroy(bam_mates_t *mates);

int bam_fetch_mate(bamFile fb, const bam_index_t *idx, int tid, int beg, 
                   int end, void *data, bam_fetch_f func, set_partition_f pfunc,
                   set_mates_f mfunc);
int samread_mate(bamFile fb, const bam_index_t *bindex, uint64_t pos0, 
                 bam_mate_iter_t *iter_p, bam_mates_t *mates);
void bam_mate_iter_destroy(bam_mate_iter_t iter);

#ifdef __cplusplus
}
#endif

#endif
