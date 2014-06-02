#ifndef PILEUPBUFFER_H
#define PILEUPBUFFER_H

#include <map>
#include <set>
#include <vector>
#include <utility>
#include <algorithm>
#include <cassert>
#include "samtools/bam.h"
#include "ResultManager.h"
#include "utilities.h"
#include <Rinternals.h>
#include <Rdefines.h>

class PileupBuffer {
  protected:
    bam_plbuf_t *plbuf;
    const char *rname;
    uint32_t start, end;
  public:
    PileupBuffer() : plbuf(NULL) {}
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
            plbuf = NULL;
        }
    }
    void plbuf_push(const bam1_t *bam) {
        bam_plbuf_push(bam, plbuf);
    }
    virtual void plbuf_init() = 0;
    virtual SEXP yield() = 0;
};

class Pileup : public PileupBuffer {
private:
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
    // cycle_bins pileupParams[9]
    int cycleBinsLength() const {
        return Rf_length(VECTOR_ELT(pileupParams, 9));
    }
    bool hasBins() const {
        return cycleBinsLength() > 0;
    }
    std::vector<int32_t> binPointsAsVec(SEXP cycleBins) const {
        int numPoints = Rf_length(cycleBins);
        std::vector<int32_t> points(numPoints);
        for(int i = 0; i != numPoints; ++i) {
            points.at(i) = INTEGER(cycleBins)[i];
        }
        return points;
    }
    int calcBin(int qpos) {
        return std::lower_bound(binPoints.begin(), binPoints.end(), qpos) -
            binPoints.begin();
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
    Pileup(SEXP schema_, SEXP pileupParams_, SEXP seqnamesLevels_)
        : schema(schema_), pileupParams(pileupParams_),
          seqnamesLevels(seqnamesLevels_), resultMgr(NULL),
          binPoints()
        {
            resultMgr =
                new ResultMgr(min_nucleotide_depth(), min_minor_allele_depth(),
                              hasStrands(), hasNucleotides(), hasBins());
            binPoints = binPointsAsVec(VECTOR_ELT(pileupParams, 9));
        }
    ~Pileup() {
        delete resultMgr;
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
        else
            Rf_error("Unrecognized nucleotide '%c'\n", nuc);
    }
};

void extract(const ResultMgrInterface * const from, SEXP to, bool hasStrands,
             bool hasNucleotides, bool hasBins);

#endif //PILEUPBUFFER_H
