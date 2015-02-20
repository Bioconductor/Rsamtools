// BamFileIterator.h:
// Iterator used when reading ranges from a bam file.

#ifndef BAMRANGEITERATOR_H
#define BAMRANGEITERATOR_H

#include "BamIterator.h"


class BamRangeIterator : public BamIterator {

    bam_iter_t iter;

    void iterate_inprogress(bamFile bfile) {
	if (iter_done)
	    return;
	if (NULL == bam) {	// first record 
	    bam = bam_init1();
	    if (bam_iter_read(bfile, iter, bam) < 0) {
		iter_done = true;
		return;
	    }
	}

	bool done = false;
	do {
	    process(bam);
	    if (bam_iter_read(bfile, iter, bam) < 0)
		iter_done = done = true;
	} while (!done);
        mate_touched_templates();
    }

    void finalize_inprogress(bamFile bfile) {
        int64_t pos = bam_tell(bfile);
        Templates::iterator it;
        // mate 'inprogress' segments for all templates
        for (it = templates.begin(); it != templates.end(); ++it)
            it->second.mate_inprogress_segments(bfile, bindex, complete,
                                                qname_prefix, qname_suffix,
                                                qname_trim);

        BamIterator::finalize_inprogress(bfile);
        bam_seek(bfile, pos, SEEK_SET);
    }

public:

    // constructor / destructor
    BamRangeIterator(const bam_index_t *bindex, int tid, int beg,
                     int end, char qname_prefix, char qname_suffix,
                     bam_qname_f qname_trim) :
        BamIterator(bindex, qname_prefix, qname_suffix, qname_trim) 
    {
	iter = bam_iter_query(bindex, tid, beg, end);
    }

    ~BamRangeIterator() {
	bam_iter_destroy(iter);
   }
};

#endif
