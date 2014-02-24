// BamFileIterator.h:
// Iterator used when reading a complete bam file.

#ifndef BAMFILEITERATOR_H
#define BAMFILEITERATOR_H

#include "BamIterator.h"


class BamFileIterator : public BamIterator {

    bool file_done;

    void iterate_inprogress(bamFile bfile) {
        if (iter_done | file_done)
            return;
        if (NULL == bam) {    // first record 
            bam = bam_init1();
            if (bam_read1(bfile, bam) < 0) {
                iter_done = true;
                return;
            }
        }

        bool done = false;
        do {
            process(bam);
            int tid = bam->core.tid;
            int pos = bam->core.pos;
            if (bam_read1(bfile, bam) < 0) {
                iter_done = file_done = done = true;
            } else if (!complete.empty()) {
                // stop if something to yield AND finished position
                done = (bam->core.tid != tid) || (bam->core.pos != pos);
            } 
        } while (!done);
        mate_touched_templates();
    }

public:

    // constructor / destructor
    BamFileIterator(const bam_index_t *bindex) :
        BamIterator(bindex), file_done(false) {}

};

#endif
