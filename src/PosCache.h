#ifndef POS_CACHE_H
#define POS_CACHE_H

#include <set>
#include <map>
#include <vector>
#include <utility>
#include <numeric>
#include <stdio.h>
#include "GenomicPosition.h"
using std::map;
using std::pair;
using std::make_pair;

struct BamTuple { // for sending info from Pileup::insert to ResultMgr
    char nuc, strand;
    int bin;
    BamTuple(char nuc_ = 'X', char strand_ = 'X', int bin_ = 0)
        : nuc(nuc_), strand(strand_), bin(bin_) { }
    bool operator<(const BamTuple& rhs) const {
        return nuc < rhs.nuc || (nuc == rhs.nuc && strand < rhs.strand) ||
            (nuc == rhs.nuc && strand == rhs.strand && bin < rhs.bin);
    }
};

struct PosCache {
    GenomicPosition genomicPosition;
    std::vector<BamTuple> tupleVec;
    std::map<char,int> nucCounts;
    typedef std::vector<BamTuple>::const_iterator tuple_iter;
    typedef std::map<char,int>::const_iterator counts_iter;
    PosCache(GenomicPosition genomicPosition_)
        : genomicPosition(genomicPosition_), tupleVec(), nucCounts() { }
    void storeTuple(BamTuple& bt) {
        tupleVec.push_back(bt);
        nucCounts.insert(make_pair(bt.nuc,0)).first->second++;
    }
    // exposed so we can use std::accumulate on nucCounts map
    static int addSecond(int i, const pair<char,int>& thePair) {
        return i + thePair.second;
    }
    int totalNucFreq() const {
        return std::accumulate(nucCounts.begin(),nucCounts.end(),0,addSecond);
    }
    int primaryNucFreq() const {
        int maxCount = 0;
        typedef map<char,int>::const_iterator iter;
        for(iter it = nucCounts.begin();
            it != nucCounts.end(); ++it) {
            if(it->second > maxCount)
                maxCount = it->second;
        }
        return maxCount;
    }
    std::set<char> passingNucs(int min) const {
        std::set<char> nucs;
        for(counts_iter it = nucCounts.begin(); it != nucCounts.end(); ++it) {
            if(it->second >= min)
                nucs.insert(it->first);
        }
        return nucs;
    }
    void clear() {
        tupleVec.clear();
        nucCounts.clear();
    }
    void print() const {
        printf("tupleVec contents:\n");
        for(tuple_iter it = tupleVec.begin(); it != tupleVec.end(); ++it) {
            printf("nuc %c str %c bin %u\n", it->nuc, it->strand, it->bin);
        }
    }
};

#endif // POS_CACHE_H
