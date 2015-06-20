#include "PileupBuffer.h"

void Pileup::signalEOI() {
    resultMgr->signalEOI();
}

int Pileup::insert(uint32_t tid, uint32_t pos, int n,
                        const bam_pileup1_t *pl, void *data)
{
    //Rprintf("pos: %d\n", pos);
    Pileup *pileup = static_cast<Pileup*>(data);
    pos = pos + 1; // 1-based indexing for R
    int bamBufOffset = 0;
    if(!pileup->isRanged || (pos >= pileup->start && pos <= pileup->end)) {
        const GenomicPosition curGenomicPosition(tid, pos);
        pileup->resultMgr->signalGenomicPosStart(curGenomicPosition);

        for(bamBufOffset = 0; bamBufOffset != n; ++bamBufOffset) {
            const bam_pileup1_t *curBam = pl + bamBufOffset;

            const uint8_t mapqual = curBam->b->core.qual;
            if(mapqual < pileup->min_mapq()) continue;

            char strand = 'X', nucleotide = 'X';
            int bin = 0;

            // positional disqualifier(s)
            if(pileup->hasBins()) { // all bin work
                const int minBinPoint = pileup->minBinPoint(),
                    maxBinPoint = pileup->maxBinPoint();
                const int qpos = curBam->qpos + 1;
                // positive bins (counting from left)
                if(minBinPoint >= 0) {
                    // discard if outside outer range
                    if(qpos > minBinPoint && qpos <= maxBinPoint)
                        bin = pileup->calcBinFromLeft(qpos);
                    else
                        continue;
                }
                // negative bins (counting from l_qseq)
                else {
                    const int32_t qlen = curBam->b->core.l_qseq;
                    // skip if qlen insufficient
                    if(qlen < (-maxBinPoint))
                        continue;
                    // distance from end
                    const int dfe = abs(qpos - qlen);
                    // use |dfe| - 1 since lowest cycle bin should be
                    // one less than leftmost pos to get the first nuc
                    if((-dfe) - 1 > minBinPoint && (-dfe) - 1 <= maxBinPoint)
                        bin = pileup->calcBinRev(dfe);
                    else
                        continue;
                }

            }

            // invariant: alignments that fail strand criterion not included
            if(pileup->hasStrands())
                strand = bam1_strand(curBam->b) ? '-' : '+';

            // IMPORTANT: it's essential that propagating insertions
            // remains in this position relative to other
            // disqualifiers; insertions are separate from called
            // bases--therefore should propagate independent of
            // considerations about the called base. I.e., an
            // insertion at a position can propagate while the called
            // base is diqualified.
            if(curBam->indel > 0 && pileup->include_insertions())
                pileup->resultMgr->forwardTuple(BamTuple('+', strand, bin));

            // whole-alignment disqualifiers/filters
            if(curBam->is_refskip) // never pileup a refskip ('N' op in cigar)
                continue;

            // individual nucleotide disqualifiers
            const uint8_t basequal = bam1_qual(curBam->b)[curBam->qpos];
            if(basequal < pileup->min_baseq()) continue;
            bool isDeletion = curBam->is_del;
            if(isDeletion && !pileup->include_deletions())
                continue;
            if(isDeletion && pileup->include_deletions())
                nucleotide = '-';
            else {
                nucleotide =
                    char(bam_nt16_rev_table[bam1_seqi(bam1_seq(curBam->b),
                                                      curBam->qpos)]);
            }

            bool dropNucleotide = (nucleotide == 'N' && pileup->ignoreNs());
            if(dropNucleotide) continue;

            pileup->resultMgr->forwardTuple(BamTuple(nucleotide, strand, bin));
        }
        const GenomicPosition maxLeftmostGenPOS(tid, pl[n-1].b->core.pos + 1);
        pileup->resultMgr->forwardLastLeftmostGenPOS(maxLeftmostGenPOS);

        // extract tuples for pos
        pileup->resultMgr->signalGenomicPosEnd();
    }
    return 0;
}

void extract(const ResultMgrInterface * const from, SEXP to, bool hasStrands,
             bool hasNucleotides, bool hasBins, bool isRanged) {
    #ifdef PILEUP_DEBUG
    assert(IS_LIST(to));
    for(int i = 0; i != Rf_length(to); ++i) {
        SEXP elt = VECTOR_ELT(to, i);
        assert(IS_LIST(elt));
        assert((unsigned int)Rf_length(elt) == from->size());
    }
    #endif
    if(!isRanged) { // seqnms must be copied from seqnmsVec
        std::copy(from->seqnmsBeg(), from->seqnmsEnd(),
                  INTEGER(VECTOR_ELT(to, 0)));
    }
    int curDim = 1; // 0 is seqnames; must start with pos
    std::copy(from->posBeg(), from->posEnd(),
              INTEGER(VECTOR_ELT(to, curDim++)));

    // Might still be worthwhile to make original vectors ints
    SEXP strand = R_NilValue, nucleotide = R_NilValue;
    if(hasStrands) {
        strand = VECTOR_ELT(to, curDim++);
        std::transform(from->strandBeg(), from->strandEnd(),
                       INTEGER(strand), Pileup::strand_to_lvl);
    }
    if(hasNucleotides) {
        nucleotide = VECTOR_ELT(to, curDim++);
        std::transform(from->nucBeg(), from->nucEnd(),
                       INTEGER(nucleotide), Pileup::nuc_to_lvl);
    }
    if(hasBins) {
        std::copy(from->binBeg(), from->binEnd(),
                  INTEGER(VECTOR_ELT(to, curDim++)));
    }

    std::copy(from->countBeg(), from->countEnd(),
              INTEGER(VECTOR_ELT(to, curDim++)));
    if(hasStrands)
        _as_strand(strand);
    if(hasNucleotides)
        _as_nucleotide(nucleotide);
}

SEXP Pileup::yield() {
    int numDims = 3;
    numDims += hasStrands() ? 1 : 0;
    numDims += hasNucleotides() ? 1 : 0;
    numDims += hasBins() ? 1 : 0;
    if(isBuffered)
        resultMgr->signalYieldStart();
    uint32_t numResults = resultMgr->size();
    SEXP result = PROTECT(Rf_allocVector(VECSXP, numDims));
    int curDim = 0;
    SET_VECTOR_ELT(result, curDim, Rf_allocVector(INTSXP, numResults));//seqns
    SEXP seqnames = VECTOR_ELT(result, curDim++);
    _as_seqlevels(seqnames, seqnamesLevels);
    // if ranged (i.e., 'which' argument passed to ScanBamParam),
    // seqnames value will be same for entire buffer otherwise, values
    // will be copied in extract function
    if(isRanged) 
        std::fill_n(INTEGER(seqnames), numResults, getSeqlevelValue(rname));
    SET_VECTOR_ELT(result, curDim++, Rf_allocVector(INTSXP, numResults)); // pos
    if(hasStrands())
        SET_VECTOR_ELT(result, curDim++, Rf_allocVector(INTSXP, numResults));
    if(hasNucleotides())
        SET_VECTOR_ELT(result, curDim++, Rf_allocVector(INTSXP, numResults));
    if(hasBins())
        SET_VECTOR_ELT(result, curDim++, Rf_allocVector(INTSXP, numResults));
    SET_VECTOR_ELT(result, curDim++, Rf_allocVector(INTSXP, numResults));//count

    SEXP nms = PROTECT(Rf_allocVector(STRSXP, numDims));
    curDim = 0;
    SET_STRING_ELT(nms, curDim++, mkChar("seqnames"));
    SET_STRING_ELT(nms, curDim++, mkChar("pos"));
    if(hasStrands())
        SET_STRING_ELT(nms, curDim++, mkChar("strand"));
    if(hasNucleotides())
        SET_STRING_ELT(nms, curDim++, mkChar("nucleotide"));
    if(hasBins())
        SET_STRING_ELT(nms, curDim++, mkChar("cycle_bin"));
    SET_STRING_ELT(nms, curDim++, mkChar("count"));
    SET_ATTR(result, R_NamesSymbol, nms);

    extract(resultMgr, result, hasStrands(), hasNucleotides(), hasBins(),
            isRanged);
    resultMgr->signalYieldEnd();

    UNPROTECT(2);
    return result;
}
