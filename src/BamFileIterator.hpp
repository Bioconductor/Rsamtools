// BamFileIterator.hpp:
// Iterator used when reading a complete bam file.

#ifndef BAMFILEITERATOR_H
#define BAMFILEITERATOR_H

#include "BamIterator.hpp"


class BamFileIterator : public BamIterator {

    uint64_t pos0;

    bool file_done;

    void iterate(bamFile bfile) {
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
            } else if (complete.size() != 0) {
                // stop if something to yield AND finished position
                done = (bam->core.tid != tid) || (bam->core.pos != pos);
            } 
        } while (!done);
    }

public:

    // constructor / destructor
    BamFileIterator(uint64_t pos0, const bam_index_t *bindex) :
        BamIterator(bindex) 
    {
        pos0 = pos0;
        file_done = false;
    }

};

#endif
