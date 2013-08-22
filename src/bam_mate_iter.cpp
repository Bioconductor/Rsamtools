#include <Rdefines.h>
#include "BamRangeIterator.hpp"
#include "BamFileIterator.hpp"
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
    mates->mates = false;
    mates->bams = NULL;
    return mates;
}

void bam_mates_realloc(bam_mates_t *mates, int n)
{
    for (int i = 0; i < mates->n; ++i) {
        bam_destroy1(mates->bams[i]);
        mates->bams[i] = NULL;
    }

    // Realloc(p, 0, *) fails inappropriately
    if (n == 0) {
	Free(mates->bams);
	mates->bams = NULL;
    } else
	mates->bams = Realloc(mates->bams, n, bam1_t *);
    mates->n = n;
    mates->mates = false;
}

void bam_mates_destroy(bam_mates_t *mates)
{
    for (int i = 0; i < mates->n; ++i)
	bam_destroy1(mates->bams[i]);
    Free(mates->bams);
    Free(mates);
}

int bam_mate_read(bamFile fb, bam_mate_iter_t iter, bam_mates_t *mates,
                  bool do_mate_all)
{
    iter->b_iter->yield(fb, mates, do_mate_all);
    return mates->n;
}

// BamRangeIterator methods
bam_mate_iter_t bam_mate_range_iter_new(const bam_index_t *bindex,
                                        int tid, int beg, int end)
{
    bam_mate_iter_t iter = Calloc(1, struct _bam_mate_iter_t);
    iter->b_iter = new BamRangeIterator(bindex, tid, beg, end);
    return iter;
}

int bam_fetch_mate(bamFile fb, const bam_index_t *idx, int tid, int beg, 
                   int end, void *data, bam_fetch_f func, set_mate_f mfunc)
{
    int i, n_rec, parse, pass;
    bam_mates_t *mates = bam_mates_new();
    bam_mate_iter_t iter = bam_mate_range_iter_new(idx, tid, beg, end);
    while ((n_rec = bam_mate_read(fb, iter, mates, true) > 0)) {
        // single yield
        pass = 0;
        for (int i = 0; i < mates->n; ++i) { 
            parse = func(mates->bams[i], data);
            if (parse == 1)
                pass = pass + 1;
        }
        if (pass > 0)
            mfunc(pass, mates->mates, data);
    }
    bam_mates_destroy(mates);
    bam_mate_iter_destroy(iter);
    return n_rec;
}

// BamFileIterator methods
bam_mate_iter_t bam_mate_file_iter_new(uint64_t pos0, const bam_index_t *bindex)
{
    bam_mate_iter_t iter = Calloc(1, struct _bam_mate_iter_t);
    iter->b_iter = new BamFileIterator(pos0, bindex);
    return iter;
}

int samread_mate(bamFile fb, const bam_index_t *bindex, uint64_t pos0,
                 bam_mate_iter_t *iter_p, bam_mates_t *mates)
{
    bam_mate_iter_t iter;
    if (NULL == *iter_p)
        *iter_p = bam_mate_file_iter_new(pos0, bindex);
    iter = *iter_p;
    iter->b_iter->iter_done = false;
    // single yield
    return bam_mate_read(fb, iter, mates, false);
}


#ifdef __cplusplus
}
#endif
