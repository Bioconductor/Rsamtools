#ifndef GENOMIC_POSITION_H
#define GENOMIC_POSITION_H

// identifier for genomic position; used for ordering of PosCache in
// associative containers
struct GenomicPosition {
    int tid, pos;
    GenomicPosition(int tid_, int pos_) : tid(tid_), pos(pos_) { }
    bool operator<(const GenomicPosition& rhs) const {
        return tid < rhs.tid || (tid == rhs.tid && pos < rhs.pos);
    }
    bool operator==(const GenomicPosition& rhs) const {
        return tid == rhs.tid && pos == rhs.pos;
    }
    bool operator<=(const GenomicPosition& rhs) const {
        return *this < rhs || *this == rhs;
    }
};

#endif // GENOMIC_POSITION_H
