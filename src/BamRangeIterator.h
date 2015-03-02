// BamRangeIterator.h:
// Iterator used when reading ranges from a bam file.

#ifndef BAMRANGEITERATOR_H
#define BAMRANGEITERATOR_H

#include "BamIterator.h"

class BamRangeIterator : public BamIterator {

    int32_t tid, beg, end;
    bam_iter_t iter;

    void iterate_inprogress(bamFile bfile) {
	if (NULL == bam) {	// first record 
	    bam = bam_init1();
	    if (bam_iter_read(bfile, iter, bam) < 0) {
		iter_done = true;
		return;
	    }
	}

	do {
	    process(bam);
	    if (bam_iter_read(bfile, iter, bam) < 0)
		iter_done = true;
	} while (!iter_done);
        mate_touched_templates();
    }

    void finalize_inprogress(bamFile bfile) {
        int64_t pos = bam_tell(bfile);
        Templates::iterator it;
        // mate 'inprogress' segments for all templates
        for (it = templates.begin(); it != templates.end(); ++it)
            it->second.mate_inprogress_segments(bfile, bindex, complete,
                                                qname_prefix, qname_suffix,
                                                tid, beg, end,
                                                header->target_len,
                                                it->first);

        BamIterator::finalize_inprogress(bfile);
        bam_seek(bfile, pos, SEEK_SET);
    }

public:

    // constructor / destructor
    BamRangeIterator(bamFile bfile, const bam_index_t *bindex,
                     int32_t tid, int32_t beg, int32_t end, 
                     BAM_DATA bam_data) :
        BamIterator(bfile, bindex, bam_data), tid(tid), beg(beg), end(end)
    {
	iter = bam_iter_query(bindex, tid, beg, end);
    }

    ~BamRangeIterator() {
	bam_iter_destroy(iter);
   }
};

#endif
