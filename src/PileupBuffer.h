#ifndef PILEUPBUFFER_H
#define PILEUPBUFFER_H

#include <map>
#include <set>
#include <vector>
#include <utility>
#include <algorithm>
#include <cassert>
#include <sam.h>
#include "ResultManager.h"
#include "utilities.h"
#include <Rinternals.h>
#include <Rdefines.h>

class PosCacheColl;

class PileupBuffer {
protected:
    bam_plbuf_t *plbuf;
    const char *rname;
    uint32_t start, end;
public:
    PileupBuffer() : plbuf((bam_plbuf_t *) NULL) {}
    virtual ~PileupBuffer() {
        plbuf_destroy();
    }
    void init(const char *_rname, const int _start, const int _end) {
        plbuf_init();
        rname = _rname;
        start = _start;
        end = _end;
    }
    void plbuf_destroy() {
        if (plbuf != NULL) {
            bam_plbuf_destroy(plbuf);
            plbuf = (bam_plbuf_t *) NULL;
        }
    }
    void plbuf_push(const bam1_t *bam) {
        bam_plbuf_push(bam, plbuf);
    }
    virtual void plbuf_init() = 0;
    virtual SEXP yield() = 0;
    virtual void signalEOI() = 0;
};

class Pileup : public PileupBuffer {
private:
    // 'ranged' means the pileup query is for specific genomic ranges;
    // 'buffered' is when intermediate results for genomic positions
    // are store across calls to pileup
    const bool isRanged, isBuffered;
    bool isQueryBin;
    int binsLength;
    const SEXP schema, pileupParams, seqnamesLevels;
    ResultMgrInterface *resultMgr;
    std::vector<int32_t> binPoints;
    int max_depth() const {
        return INTEGER(VECTOR_ELT(pileupParams, 0))[0];
    }
    uint8_t min_baseq() const {
        return INTEGER(VECTOR_ELT(pileupParams, 1))[0];
    }
    uint8_t min_mapq() const {
        return INTEGER(VECTOR_ELT(pileupParams, 2))[0];
    }
    int min_nucleotide_depth() const {
        return INTEGER(VECTOR_ELT(pileupParams, 3))[0];
    }
    int min_minor_allele_depth() const {
        return INTEGER(VECTOR_ELT(pileupParams, 4))[0];
    }
    bool hasStrands() const {
        return LOGICAL(VECTOR_ELT(pileupParams, 5))[0]; // distinguish_strands
    }
    bool hasNucleotides() const {
        return LOGICAL(VECTOR_ELT(pileupParams, 6))[0]; // distinguish_nucs
    }
    bool ignoreNs() const {
        return LOGICAL(VECTOR_ELT(pileupParams, 7))[0];
    }
    bool include_deletions() const {
        return LOGICAL(VECTOR_ELT(pileupParams, 8))[0];
    }
    bool include_insertions() const {
        return LOGICAL(VECTOR_ELT(pileupParams, 9))[0];
    }
    // left_bins pileupParams[10]
    int leftBinsParamLength() const {
        return Rf_length(VECTOR_ELT(pileupParams, 10));
    }
    // query_bins pileupParams[11]
    int queryBinsParamLength() const {
        return Rf_length(VECTOR_ELT(pileupParams, 11));
    }
    int getBinsLength() const {
        return binsLength;
    }
    bool hasBins() const {
        return binsLength > 0;
    }
    bool isQueryBinMode() const {
        return isQueryBin;
    }
    std::vector<int32_t> binPointsAsVec(SEXP cycleBins) const {
        int numPoints = Rf_length(cycleBins);
        std::vector<int32_t> points(numPoints);
        for(int i = 0; i != numPoints; ++i) {
            points.at(i) = INTEGER(cycleBins)[i];
        }
        return points;
    }
    int calcBin(int dfe) {
        return std::lower_bound(binPoints.begin(), binPoints.end(), dfe)
            - binPoints.begin();
    }
    int32_t maxBinPoint() const {
        return *(binPoints.end() - 1);
    }
    int32_t minBinPoint() const {
        return *binPoints.begin();
    }
    int getSeqlevelValue(const char* theRname) const {
        int idx = 0;
        for(idx = 0; idx != Rf_length(seqnamesLevels); ++idx) {
            const char* curLevel = CHAR(STRING_ELT(seqnamesLevels, idx));
            //Rprintf("rname (%s) =? level (%s)\n", rname, curLevel);
            if(strcmp(theRname, curLevel) == 0)
                return idx + 1;
        }
        Rf_error("rname '%s' not in seqnames levels", rname);
        return -1;
    }
public:
    // void *posCacheColl is dumbly passed through to ResultMgr to
    // give ResultMgr pointer to the struct _BAM_FILE pbuffer member;
    // PileupBuffer doesn't know anything about PosCacheColl
    Pileup(bool isRanged_, bool isBuffered_, SEXP schema_, SEXP pileupParams_,
           SEXP seqnamesLevels_, PosCacheColl** posCacheColl_)
        : isRanged(isRanged_), isBuffered(isBuffered_), isQueryBin(false),
          binsLength(0),
          schema(schema_),
          pileupParams(pileupParams_), seqnamesLevels(seqnamesLevels_),
          resultMgr((ResultMgrInterface *) NULL), binPoints()
        {
            if(isRanged && isBuffered) {
                Rf_error("internal: Pileup cannot both query specific genomic ranges and store partial genomic position results");
            }
            // bins stuff (needs to come before ResultMgr ctor call
            // because ctor uses hasBins()
            if(leftBinsParamLength() > 0) {
                isQueryBin = false;
                binsLength = leftBinsParamLength();
                binPoints = binPointsAsVec(VECTOR_ELT(pileupParams, 10));
            } else if(queryBinsParamLength() > 0) {
                isQueryBin = true;
                binsLength = queryBinsParamLength();
                binPoints = binPointsAsVec(VECTOR_ELT(pileupParams, 11));
            }
            resultMgr =
                new ResultMgr(min_nucleotide_depth(), min_minor_allele_depth(),
                              hasStrands(), hasNucleotides(), hasBins(),
                              isRanged, isBuffered, posCacheColl_);
        }
    ~Pileup() {
        delete resultMgr;
    }
    bool needMoreInput() const {
        bool needMore = resultMgr->numYieldablePosCaches() == 0;
        return needMore;
    }
    bool isBufferedPileup() const {
        return isBuffered;
    }
    void plbuf_init() {
        plbuf = bam_plbuf_init(insert, this);
        int theDepth = max_depth();
        if(theDepth < 1)
            Rf_error("'max_depth' must be greater than 0, got '%d'", theDepth);
        // +1 essential because when maxcnt > 1 num reads processed = maxcnt - 1
        int num_reads_to_process = theDepth < 2 ? 1 : theDepth + 1;
        bam_plp_set_maxcnt(plbuf->iter, num_reads_to_process);
    }
    static int insert(uint32_t tid, uint32_t pos, int n,
                      const bam_pileup1_t *pl, void *data);
    SEXP yield();
    void signalEOI();
    static int strand_to_lvl(char strand) {
        return strand == '+' ? 1 : 2;
    }
    static int nuc_to_lvl(char nuc) {
        if(nuc == 'A')
            return 1;
        else if(nuc == 'C')
            return 2;
        else if(nuc == 'G')
            return 3;
        else if(nuc == 'T')
            return 4;
        else if(nuc == 'N')
            return 5;
        else if(nuc == '=')
            return 6;
        else if(nuc == '-')
            return 7;
        else if(nuc == '+')
            return 8;
        else
            Rf_error("Unrecognized nucleotide '%c'\n", nuc);
    }
};

void extract(const ResultMgrInterface * const from, SEXP to, bool hasStrands,
             bool hasNucleotides, bool hasBins);

#endif //PILEUPBUFFER_H
