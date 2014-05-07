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
        /* FIXME: separate initialize method when space == NULL */
        if (R_NilValue == space) {
            Rf_error("missing 'which' not yet implemented");
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
    void plbuf_push(const bam1_t *bam) {
        buffer.plbuf_push(bam);
    }
};

#endif // PILEUPBUFFERSHIM_H
