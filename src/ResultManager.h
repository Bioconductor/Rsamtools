#ifndef RESULT_MANAGER_H
#define RESULT_MANAGER_H

#include <map>
#include <vector>
#include <utility>
#ifdef PILEUP_DEBUG
#include <stdio.h>
#include <cassert>
#include "nate_utilities.h"
#endif
#include "PosCache.h"
#include "PosCacheColl.h"
#include "GenomicPosition.h"
#include <Rinternals.h>
typedef std::vector<int>::const_iterator int_const_it;
typedef std::vector<char>::const_iterator char_const_it;

class ResultMgrInterface
{
public:
    virtual void signalGenomicPosStart(const GenomicPosition& genPos) = 0;
    virtual void forwardLastLeftmostGenPOS(const GenomicPosition& genPos) = 0;
    virtual void forwardTuple(BamTuple bTuple) = 0;
    virtual void extractFromPosCache() = 0;
    virtual void signalGenomicPosEnd() = 0;
    virtual int size() const = 0;
    virtual void signalYieldStart() = 0;
    virtual void signalYieldEnd() = 0;
    virtual int numYieldablePosCaches() const = 0;
    virtual void signalEOI() = 0;
    virtual ~ResultMgrInterface() {}
    virtual int_const_it seqnmsBeg() const = 0;
    virtual int_const_it seqnmsEnd() const = 0;
    virtual int_const_it posBeg() const = 0;
    virtual int_const_it posEnd() const = 0;
    virtual int_const_it countBeg() const = 0;
    virtual int_const_it countEnd() const = 0;
    virtual char_const_it strandBeg() const = 0;
    virtual char_const_it strandEnd() const = 0;
    virtual char_const_it nucBeg() const = 0;
    virtual char_const_it nucEnd() const = 0;
    virtual int_const_it binBeg() const = 0;
    virtual int_const_it binEnd() const = 0;
};

class ResultMgr : public ResultMgrInterface
{
private:
    std::vector<int> seqnmsVec, posVec, binVec, countVec;
    std::vector<char> strandVec, nucVec;
    PosCache* posCache;
    // posCacheCollptrptr is ptr to struct _BAM_FILE's pbuffer ptr
    PosCacheColl** posCacheCollptrptr;
    const int min_nuc_depth, min_minor_allele_depth;
    const bool hasStrands, hasNucleotides, hasBins, isRanged, isBuffered;
    // lastLeftmostGenPOS is for bookkeeping for buffered pileups;
    // PosCaches that correspond to completed positions for which
    // posCache.genomicPosition < minLeftmostGenPOS are completed and,
    // therefore, ready to be extracted
    GenomicPosition lastLeftmostGenPOS; // tid + POS
    // posCachePassesFilters:
    // postondition: returns true if posCache's object satisfies the
    // filtering criteria that are applied to completed positions
    bool posCachePassesFilters(const PosCache& posCachePtr);
    template <bool wantNuc, bool wantStrand, bool wantBin>
    void doExtractFromPosCache();
public:
    // used only in buffered and unbuffered whole file contexts??
    virtual void signalGenomicPosStart(const GenomicPosition& genPos);
    void forwardLastLeftmostGenPOS(const GenomicPosition& genPos);
    virtual void forwardTuple(BamTuple bTuple);
    virtual void extractFromPosCache();
    template <bool wantNuc, bool wantStrand, bool wantBin>
    void doExtractFromPosCache(const std::set<char>& nucs);
    virtual void signalGenomicPosEnd();
    virtual int size() const;
    virtual void printVecs() const;
    void signalYieldStart();
    virtual void signalYieldEnd();
    virtual int numYieldablePosCaches() const;
    virtual void signalEOI();
    virtual ~ResultMgr() {
        // FIX ME: must deallocate posCacheColl only if done with it!!!
        //delete posCacheColl;
        //posCacheColl = NULL;
    }
    ResultMgr(int min_nuc_depth_, int min_minor_allele_depth_,
              bool hasStrands_, bool hasNucleotides_, bool hasBins_,
              bool isRanged_, bool isBuffered_, PosCacheColl** posCacheColl_) :
        seqnmsVec(), posVec(),  binVec(), countVec(), strandVec(), nucVec(),
        posCache(), posCacheCollptrptr(posCacheColl_),
        min_nuc_depth(min_nuc_depth_),
        min_minor_allele_depth(min_minor_allele_depth_),
        hasStrands(hasStrands_), hasNucleotides(hasNucleotides_),
        hasBins(hasBins_), isRanged(isRanged_), isBuffered(isBuffered_),
        lastLeftmostGenPOS(0, 0)
        {
            if(isBuffered && *posCacheCollptrptr == NULL) {
                *posCacheCollptrptr = new PosCacheColl();
            }
        }
    virtual int_const_it seqnmsBeg() const;
    virtual int_const_it seqnmsEnd() const;
    virtual int_const_it posBeg() const;
    virtual int_const_it posEnd() const;
    virtual int_const_it binBeg() const;
    virtual int_const_it binEnd() const;
    virtual int_const_it countBeg() const;
    virtual int_const_it countEnd() const;
    virtual char_const_it strandBeg() const;
    virtual char_const_it strandEnd() const;
    virtual char_const_it nucBeg() const;
    virtual char_const_it nucEnd() const;
};

#endif // RESULT_MANAGER_H
