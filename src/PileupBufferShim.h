#ifndef PILEUPBUFFERSHIM_H
#define PILEUPBUFFERSHIM_H

#include "PileupBuffer.h"

class PileupBufferShim {

private:
    const SEXP space;
    SEXP result;
    PileupBuffer &buffer;
public:
    PileupBufferShim(SEXP _space, SEXP _result, PileupBuffer &_buffer) :
        space(_space), result(_result), buffer(_buffer)
        {}

    void start1(const int irange) {
        if (R_NilValue == space) {
            buffer.init(NULL, 0, 0);
        } else {
            buffer.init(
                CHAR(STRING_ELT(VECTOR_ELT(space, 0), irange)),
                INTEGER(VECTOR_ELT(space, 1))[irange],
                INTEGER(VECTOR_ELT(space, 2))[irange]);
        }
    }
    void finish1(const int irange) {
        plbuf_push(0);
        SET_VECTOR_ELT(result, irange, buffer.yield());
        buffer.plbuf_destroy();
    }
    // The only way to trigger running the callback function
    // (Pileup::insert in this case) is to push NULL to the buffer and
    // destroy it. Therefore, must destroy and recreate plbuf each
    // time yieldSize records are pushed.
    void process_yieldSize_chunk() {
        plbuf_push(0);
        buffer.plbuf_destroy(); // trigger run of Pileup::insert
        buffer.init(NULL, 0, 0);
    }
    // intended to be called from _pileup_bam after EOI message sent
    // to PileupBuffer; same as finish1 but no plbuf is in use
    void flush() {
        //Rprintf("flushing\n");
        SET_VECTOR_ELT(result, 0, buffer.yield());
    }
    void plbuf_push(const bam1_t *bam) {
        buffer.plbuf_push(bam);
    }
};

#endif // PILEUPBUFFERSHIM_H
