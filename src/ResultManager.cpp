#include "ResultManager.h"

void ResultMgr::signalPosStart(int pos) {
    curPos = pos;
}

void ResultMgr::forwardTuple(char strand, char nucleotide) {
    posCache.increment(strand, nucleotide);
}

void ResultMgr::signalPosEnd() {
    typedef map<char,StrandTree>::const_iterator nuc_iter;
    typedef map<char,int>::const_iterator strand_iter;
    //printf("\nsignalPosEnd:\n");
    int totalNucCount = posCache.count();
    //printf("total nuc pop: %d\n", totalNucCount);
    int primaryNucCount = posCache.primaryNucCount();
    //printf("primaryNucCount %d\n", primaryNucCount);
    // must always check min_minor_allele_depth
    if(totalNucCount - primaryNucCount >= min_minor_allele_depth) {
        if(hasStrands && !hasNucleotides)
            extractStrandOnly();
        else if(!hasStrands && !hasNucleotides)
            extractCollapseAll();
        else {
            for(nuc_iter nuc_it = posCache.nucMap.begin();
                nuc_it != posCache.nucMap.end(); ++nuc_it) {
                char nuc = nuc_it->first;
                //printf("nuc %c nucCount %u\n", nuc, nuc_it->second.count());
                if(nuc_it->second.count() < min_nuc_depth) {
                    //printf("Skipping nucleotide\n");
                    continue;
                }
                for(strand_iter strand_it = nuc_it->second.strandMap.begin();
                    strand_it != nuc_it->second.strandMap.end(); ++strand_it) {
                    char strand = strand_it->first;
                    int count = strand_it->second;
                    //printf("  pos %u nuc %c strand %c cnt %u\n",curPos,nuc,strand,count);
                    posVec.push_back(curPos);
                    countVec.push_back(count);
                    if(hasStrands)
                        strandVec.push_back(strand);
                    if(hasNucleotides)
                        nucVec.push_back(nuc);
                }
            }
        }
    }
    posCache.clear();

    // debugging
    // nate_uts::printcoll(posVec);
    // nate_uts::printcoll(countVec);
    // if(hasStrands) nate_uts::printcoll(strandVec);
    // if(hasNucleotides) nate_uts::printcoll(nucVec);
}

void ResultMgr::extractStrandOnly() {
    typedef map<char,StrandTree>::const_iterator nuc_iter;
    typedef map<char,int>::const_iterator strand_iter;
    int plusCount = 0, minusCount = 0;
    for(nuc_iter nuc_it = posCache.nucMap.begin();
        nuc_it != posCache.nucMap.end(); ++nuc_it) {
        for(strand_iter strand_it = nuc_it->second.strandMap.begin();
            strand_it != nuc_it->second.strandMap.end(); ++strand_it) {
            char strand = strand_it->first;
            assert(strand == '+' || strand == '-');
            if(strand == '+') plusCount += strand_it->second;
            else minusCount += strand_it->second;
        }
    }
    if(plusCount > 0) {
        posVec.push_back(curPos);
        strandVec.push_back('+');
        countVec.push_back(plusCount);
    }
    if(minusCount > 0) {
        posVec.push_back(curPos);
        strandVec.push_back('-');
        countVec.push_back(minusCount);
    }
}

// this special case required (instead of just using posCache.count())
// because need to enforce min_nuc_depth criterion on each nucleotide
void ResultMgr::extractCollapseAll() {
    typedef map<char,StrandTree>::const_iterator nuc_iter;
    int posTotal = 0; // across all nucs that pass min_nuc_depth criterion
    for(nuc_iter nuc_it = posCache.nucMap.begin();
        nuc_it != posCache.nucMap.end(); ++nuc_it) {
        // don't care which nuc
        if(nuc_it->second.count() < min_nuc_depth) {
            //printf("Skipping nucleotide\n");
            continue;
        }
        posTotal += nuc_it->second.count();
    }
    if(posTotal > 0) {
        posVec.push_back(curPos);
        countVec.push_back(posTotal);
    }
}

int ResultMgr::size() const {
    // number of tuples in linear representation
    return posVec.size();
}

void ResultMgr::signalYieldEnd() {
    posVec.clear();
    countVec.clear();
    strandVec.clear();
    nucVec.clear();
}

int_const_it ResultMgr::posBeg() const { return posVec.begin(); }
int_const_it ResultMgr::posEnd() const { return posVec.end(); }
int_const_it ResultMgr::countBeg() const { return countVec.begin(); }
int_const_it ResultMgr::countEnd() const { return countVec.end(); }
char_const_it ResultMgr::strandBeg() const { return strandVec.begin(); }
char_const_it ResultMgr::strandEnd() const { return strandVec.end(); }
char_const_it ResultMgr::nucBeg() const { return nucVec.begin(); }
char_const_it ResultMgr::nucEnd() const { return nucVec.end(); }
