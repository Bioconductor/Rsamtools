#include <Rdefines.h>
#include "BamRangeIterator.h"
#include "BamFileIterator.h"
#include "bam_mate_iter.h"

#ifdef __cplusplus
extern "C" {
#endif

struct _bam_mate_iter_t {
    BamIterator *b_iter;
};

// BamIterator methods
void bam_mate_iter_destroy(bam_mate_iter_t iter)
{
    delete iter->b_iter;
    Free(iter);
}

bam_mates_t *bam_mates_new()
{
    bam_mates_t *mates = Calloc(1, bam_mates_t);
    mates->n = 0;
    mates->mated = MATE_UNKNOWN;
    mates->bams = (const bam1_t **) NULL;
    return mates;
}

void bam_mates_realloc(bam_mates_t *result, int n, MATE_STATUS mated)
{
    for (int i = 0; i < result->n; ++i) {
        bam_destroy1((bam1_t *) result->bams[i]);
        result->bams[i] = (const bam1_t *) NULL;
    }

    // Realloc(p, 0, *) fails inappropriately
    if (n == 0) {
        Free(result->bams);
        result->bams = (const bam1_t **) NULL;
    } else {
        result->bams = Realloc(result->bams, n, const bam1_t *);
    }
    result->n = n;
    result->mated = mated;
}

void bam_mates_destroy(bam_mates_t *mates)
{
    for (int i = 0; i < mates->n; ++i)
        bam_destroy1((bam1_t *) mates->bams[i]);
    Free(mates->bams);
    Free(mates);
}

int bam_mate_read(bamFile fb, bam_mate_iter_t iter, bam_mates_t *mates)
{
    iter->b_iter->yield(fb, mates);
    return mates->n;
}

// BamRangeIterator methods
bam_mate_iter_t bam_mate_range_iter_new(bamFile bfile,
                                        const bam_index_t *bindex, int tid,
                                        int beg, int end)
{
    bam_mate_iter_t iter = Calloc(1, struct _bam_mate_iter_t);
    iter->b_iter = new BamRangeIterator(bfile, bindex, tid, beg, end);
    return iter;
}

int bam_fetch_mate(bamFile bf, const bam_index_t *idx, int tid, int beg, 
                   int end, void *data, bam_fetch_mate_f func)
{
    BAM_DATA bd = (BAM_DATA) data;
    int n_rec;
    bam_mates_t *mates = bam_mates_new();
    bam_mate_iter_t iter = bam_mate_range_iter_new(bf, idx, tid, beg, end);
    iter->b_iter->set_bam_data(bd);
    while ((n_rec = bam_mate_read(bf, iter, mates) > 0))
        func(mates, data);
    bam_mate_iter_destroy(iter);
    bam_mates_destroy(mates);
    return n_rec;
}

// BamFileIterator methods
bam_mate_iter_t bam_mate_file_iter_new(bamFile bfile,
                                       const bam_index_t *bindex)
{
    bam_mate_iter_t iter = Calloc(1, struct _bam_mate_iter_t);
    iter->b_iter = new BamFileIterator(bfile, bindex);
    return iter;
}

int samread_mate(bamFile bfile, const bam_index_t *bindex,
                 bam_mate_iter_t *iter_p, bam_mates_t *mates,
                 void *data)
{
    BAM_DATA bd = (BAM_DATA) data;
    bam_mate_iter_t iter;
    int status;
    if ((bam_mate_iter_t) NULL == *iter_p)
        *iter_p = bam_mate_file_iter_new(bfile, bindex);
    iter = *iter_p;
    iter->b_iter->set_bam_data(bd);
    iter->b_iter->iter_done = false;
    // single yield
    status = bam_mate_read(bfile, iter, mates);
    iter->b_iter->set_bam_data((BAM_DATA) NULL);
    return status;
}


#ifdef __cplusplus
}
#endif
