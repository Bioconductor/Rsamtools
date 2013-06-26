#include <Rdefines.h>
#include "BamRangeIterator.hpp"
#include "bam_mate_iter.h"

#ifdef __cplusplus
extern "C" {
#endif

struct _bam_mate_iter_t {
    BamRangeIterator *b_iter;
};

bam_mate_iter_t bam_mate_query(const bam_index_t *bindex,
                               int tid, int beg, int end)
{
    bam_mate_iter_t iter = Calloc(1, struct _bam_mate_iter_t);
    iter->b_iter = new BamRangeIterator(bindex, tid, beg, end, true);
    return iter;
}

int bam_mate_read(bamFile fp, bam_mate_iter_t iter, bam_mates_t *mates)
{
    list<bam1_t *> elts = iter->b_iter->yield(fp);
    bam_mates_realloc(mates, elts.size());
    int i = 0;
    while (!elts.empty()) {
	mates->bams[i++] = elts.front();
	elts.pop_front();
    }
    return mates->n != 0;
}

void bam_mate_destroy(bam_mate_iter_t iter)
{
    delete iter->b_iter;
    Free(iter);
}

bam_mates_t *bam_mates_new()
{
    bam_mates_t *mates = Calloc(1, bam_mates_t);
    mates->n = 0;
    mates->bams = NULL;
    return mates;
}

void bam_mates_realloc(bam_mates_t *mates, int n)
{
    for (int i = 0; i < mates->n; ++i)
        bam_destroy1(mates->bams[i]);

    // Realloc(p, 0, *) fails inappropriately
    if (n == 0) {
	Free(mates->bams);
	mates->bams = NULL;
    } else
	mates->bams = Realloc(mates->bams, n, bam1_t *);
    mates->n = n;
}

void bam_mates_destroy(bam_mates_t *mates)
{
    for (int i = 0; i < mates->n; ++i)
	bam_destroy1(mates->bams[i]);
    Free(mates->bams);
    Free(mates);
}

int bam_mate_fetch(bamFile fb, const bam_index_t *idx, int tid, int beg, 
                   int end, void *data, bam_fetch_f func)
{
    int n_rec;
    bam_mates_t *result = bam_mates_new();
    bam_mate_iter_t iter = bam_mate_query(idx, tid, beg, end);
    while (n_rec = bam_mate_read(fb, iter, result) > 0)
        for (int i = 0; i < result->n; ++i)
            func(result->bams[i], data);
    bam_mate_destroy(iter);
    bam_mates_destroy(result);
    return n_rec;
}

int samread_mate(samfile_t *fp, bam_mates_t *b)
{
    return 0;
}

#ifdef __cplusplus
}
#endif
