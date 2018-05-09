#ifndef BAM_MATE_ITER_H
#define BAM_MATE_ITER_H

#include <sam.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    MATE_UNKNOWN, MATE_MATED, MATE_AMBIGUOUS, MATE_UNMATED
} MATE_STATUS;

typedef struct _bam_mate_iter_t *bam_mate_iter_t;

typedef struct {
    const bam1_t **bams;
    int n;
    MATE_STATUS mated;
} bam_mates_t;

typedef int (*bam_fetch_mate_f)(const bam_mates_t *b, void *data);

bam_mates_t *bam_mates_new();
void bam_mates_realloc(bam_mates_t *mates, int n, MATE_STATUS mated);
void bam_mates_destroy(bam_mates_t *mates);

int bam_fetch_mate(bamFile fb, const bam_index_t *idx, int tid, int beg, 
                   int end, void *data, bam_fetch_mate_f func);
int samread_mate(bamFile fb, const bam_index_t *bindex,
                 bam_mate_iter_t *iter_p, bam_mates_t *mates,
                 void *data);
void bam_mate_iter_destroy(bam_mate_iter_t iter);

#ifdef __cplusplus
}
#endif

#endif
