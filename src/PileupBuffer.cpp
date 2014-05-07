#include "PileupBuffer.h"

int Pileup::insert(uint32_t tid, uint32_t pos, int n,
                        const bam_pileup1_t *pl, void *data)
{
    //Rprintf("pos: %d\n", pos);
    Pileup *pileup = (Pileup *) data;
    pos = pos + 1;
    int bamBufOffset = 0;
    if(pos >= pileup->start && pos <= pileup->end) {
        pileup->resultMgr->signalPosStart(pos);
        for(bamBufOffset = 0; bamBufOffset != n; ++bamBufOffset) {
            const bam_pileup1_t *curBam = pl + bamBufOffset;
            // int nucCode = bam1_seqi(bam1_seq(curBam->b), curBam->qpos);
            // Rprintf("pos %d nuc code %d\n", pos, nucCode);
            // Rprintf("\tchar %c\n", bam_nt16_rev_table[nucCode]);
            

            // whole-alignment disqualifiers/filters
            if(curBam->is_refskip) // never pileup a refskip ('N' op in cigar)
                continue;
            const uint8_t mapqual = curBam->b->core.qual;
            if(mapqual < pileup->min_mapq()) continue;
            
            // individual nucleotide disqualifiers
            char strand = 'X', nucleotide = 'X';
            const uint8_t basequal = bam1_qual(curBam->b)[curBam->qpos];
            if(basequal < pileup->min_baseq()) continue;
            // invariant: alignments that fail strand criterion not included
            if(pileup->hasStrands())
                strand = bam1_strand(curBam->b) ? '-' : '+';
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
            // remove following line for production
            //int nucCode = bam1_seqi(bam1_seq(curBam->b), curBam->qpos);
            //Rprintf("strand: %c nuc code: %d nuc: %c\n", strand, nucCode, nucleotide);
            bool dropNucleotide = (nucleotide == 'N' && pileup->ignoreNs());
            if(dropNucleotide) continue;

            pileup->resultMgr->forwardTuple(strand, nucleotide);
        }
        // extract tuples for pos
        pileup->resultMgr->signalPosEnd();
    }
    return 0;
}

void extract(const ResultMgrInterface * const from, SEXP to, bool hasStrands,
    bool hasNucleotides) {
    #ifdef PILEUP_DEBUG
    assert(IS_LIST(to));
    for(int i = 0; i != Rf_length(to); ++i) {
        SEXP elt = VECTOR_ELT(to, i);
        assert(IS_LIST(elt));
        assert((unsigned int)Rf_length(elt) == from->size());
    }
    #endif
    int curDim = 1; // 0 is seqnames; must start with pos
    std::copy(from->posBeg(), from->posEnd(),
              INTEGER(VECTOR_ELT(to, curDim++)));

    SEXP strand = R_NilValue, nucleotide = R_NilValue;
    if(hasStrands)
        strand = VECTOR_ELT(to, curDim++);
    if(hasNucleotides)
        nucleotide = VECTOR_ELT(to, curDim++);
    char_const_it strandIt = from->strandBeg(), nucIt = from->nucBeg();
    // Might be better to use integers in original vector so can use
    // std::copy algorithm; it's like this now because of storing the
    // chars; alternatively, use std::transform
    for(int i = 0; strandIt != from->strandEnd() || nucIt != from->nucEnd();
        ++i) {
        // *pos++ = posVec.at(i);
        // *count++ = countVec.at(i);
        if(hasStrands) {
            // const char *theStrand = &(*strandIt);
            // SET_STRING_ELT(strand, i, mkCharLen(theStrand, 1));
            int theStrand = *strandIt == '+' ? 1 : 2;
            INTEGER(strand)[i] = theStrand;
            ++strandIt;
        }
        if(hasNucleotides) {
            int theNuc = Pileup::nuc_to_lvl(*nucIt);
            //Rprintf("theNuc value: %d\n", theNuc);
            INTEGER(nucleotide)[i] = theNuc;
            ++nucIt;
        }
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
    uint32_t numResults = resultMgr->size();
    SEXP result = PROTECT(Rf_allocVector(VECSXP, numDims));
    int curDim = 0;
    SET_VECTOR_ELT(result, curDim, Rf_allocVector(INTSXP, numResults));//seqns
    SEXP seqnames = VECTOR_ELT(result, curDim++);
    std::fill_n(INTEGER(seqnames), numResults, getSeqlevelValue(rname));
    SET_VECTOR_ELT(result, curDim++, Rf_allocVector(INTSXP, numResults)); // pos
    if(hasStrands())
        SET_VECTOR_ELT(result, curDim++, Rf_allocVector(INTSXP, numResults));
    if(hasNucleotides())
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
    SET_STRING_ELT(nms, curDim++, mkChar("count"));
    SET_ATTR(result, R_NamesSymbol, nms);

    extract(resultMgr, result, hasStrands(), hasNucleotides());
    resultMgr->signalYieldEnd();
    _as_seqlevels(VECTOR_ELT(result, 0), seqnamesLevels);

    UNPROTECT(2);
    return result;
}
