#include "ResultManager.h"

void ResultMgr::signalGenomicPosStart(const GenomicPosition& genPos) {
    //Rprintf("signalGenomicPosStart tid %d pos %d\n", genPos.tid, genPos.pos);
    if(isBuffered && (posCache != (PosCache *) NULL))
        Rf_error("internal: ResultMgr's previous posCache not deallocated");
    posCache = new PosCache(genPos);
    if(isBuffered) {
        getPosCacheFromColl(**posCacheCollptrptr, posCache);
    }
    //posCache->print();
    //Rprintf("end of signalGenomicPosStart\n\n");
}

void ResultMgr::forwardLastLeftmostGenPOS(const GenomicPosition& genPos) {
    //Rprintf("forwardMinLeftmostGenPos tid %d pos %d\n", genPos.tid, genPos.pos);
    lastLeftmostGenPOS = genPos;
}

void ResultMgr::forwardTuple(BamTuple bTuple) {
    posCache->storeTuple(bTuple);
}

template <bool wantNuc, bool wantStrand, bool wantBin>
void ResultMgr::doExtractFromPosCache(const std::set<char>&) { }

template <> void // distinguish nothing / collapse all
ResultMgr::doExtractFromPosCache<false,false,false>(const std::set<char>& nucs){
    typedef std::map<char,int>::const_iterator iter;
    int count = 0;
    for(iter it = posCache->nucCounts.begin(); it != posCache->nucCounts.end();
        ++it) {
        if(nucs.find(it->first) != nucs.end()) {
            count += it->second;
        }
    }
    if(count > 0) {
        countVec.push_back(count);
    }
}

template <> void // maximally distinguished
ResultMgr::doExtractFromPosCache<true,true,true>(const std::set<char>& nucs) {
    std::map<BamTuple,int> tupleMap;
    typedef std::vector<BamTuple>::const_iterator iter;
    for(iter it = posCache->tupleVec.begin(); it != posCache->tupleVec.end();
        ++it) {
        if(nucs.find(it->nuc) != nucs.end())
            tupleMap.insert(make_pair(*it, 0)).first->second++;
    }
    typedef std::map<BamTuple,int>::const_iterator count_iter;
    for(count_iter it = tupleMap.begin(); it != tupleMap.end(); ++it) {
        countVec.push_back(it->second);
        nucVec.push_back(it->first.nuc);
        strandVec.push_back(it->first.strand);
        binVec.push_back(it->first.bin);
    }
}

template <> void // A (nuc only)
ResultMgr::doExtractFromPosCache<true,false,false>(const std::set<char>& nucs){
    typedef std::map<char,int>::const_iterator iter;
    for(iter it = posCache->nucCounts.begin(); it != posCache->nucCounts.end();
        ++it) {
        if(nucs.find(it->first) != nucs.end()) {
            countVec.push_back(it->second);
            nucVec.push_back(it->first);
        }
    }
}

template <> void // B (strand only)
ResultMgr::doExtractFromPosCache<false,true,false>(const std::set<char>& nucs){
    std::map<char,int> strandMap; //<strand,count>
    typedef std::vector<BamTuple>::const_iterator iter;
    for(iter it = posCache->tupleVec.begin(); it != posCache->tupleVec.end();
        ++it) {
        if(nucs.find(it->nuc) != nucs.end()) {
            strandMap.insert(make_pair(it->strand,0)).first->second++;
        }
    }
    typedef std::map<char,int>::const_iterator char_iter;
    for(char_iter it = strandMap.begin(); it != strandMap.end(); ++it) {
        countVec.push_back(it->second);
        strandVec.push_back(it->first);
    }
}

template <> void // C (bin only)
ResultMgr::doExtractFromPosCache<false,false,true>(const std::set<char>& nucs){
    std::map<int,int> binMap; //<bin,count>
    typedef std::vector<BamTuple>::const_iterator iter;
    for(iter it = posCache->tupleVec.begin(); it != posCache->tupleVec.end();
        ++it) {
        if(nucs.find(it->nuc) != nucs.end()) {
            binMap.insert(make_pair(it->bin,0)).first->second++;
        }
    }
    typedef std::map<int,int>::const_iterator count_iter;
    for(count_iter it = binMap.begin(); it != binMap.end(); ++it) {
        countVec.push_back(it->second);
        binVec.push_back(it->first);
    }
}

template <> void // AB (nuc and strand)
ResultMgr::doExtractFromPosCache<true,true,false>(const std::set<char>& nucs){
    std::map<pair<char,char>,int> nucstrMap; // <<nuc,strand>,count>
    typedef std::vector<BamTuple>::const_iterator iter;
    for(iter it = posCache->tupleVec.begin(); it != posCache->tupleVec.end();
        ++it) {
        if(nucs.find(it->nuc) != nucs.end()) {
            pair<char,char> tmp = make_pair(it->nuc,it->strand);
            nucstrMap.insert(make_pair(tmp,0)).first->second++;
        }
    }
    typedef std::map<pair<char,char>,int>::const_iterator count_iter;
    for(count_iter it = nucstrMap.begin(); it != nucstrMap.end(); ++it) {
        countVec.push_back(it->second);
        nucVec.push_back(it->first.first);
        strandVec.push_back(it->first.second);
    }
}

template <> void // BC (strand and bin)
ResultMgr::doExtractFromPosCache<false,true,true>(const std::set<char>& nucs){
    std::map<pair<char,int>,int> strbinMap; // <<strand,bin>,count>
    typedef std::vector<BamTuple>::const_iterator iter;
    for(iter it = posCache->tupleVec.begin(); it != posCache->tupleVec.end();
        ++it) {
        if(nucs.find(it->nuc) != nucs.end()) {
            pair<char,int> tmp = make_pair(it->strand,it->bin);
            strbinMap.insert(make_pair(tmp,0)).first->second++;
        }
    }
    typedef std::map<pair<char,int>,int>::const_iterator count_iter;
    for(count_iter it = strbinMap.begin(); it != strbinMap.end(); ++it) {
        countVec.push_back(it->second);
        strandVec.push_back(it->first.first);
        binVec.push_back(it->first.second);
    }
}

template <> void // AC (nuc and bin)
ResultMgr::doExtractFromPosCache<true,false,true>(const std::set<char>& nucs){
    std::map<pair<char,int>,int> nucbinMap; // <<nuc,bin>,count>
    typedef std::vector<BamTuple>::const_iterator iter;
    for(iter it = posCache->tupleVec.begin(); it != posCache->tupleVec.end();
        ++it) {
        if(nucs.find(it->nuc) != nucs.end()) {
            pair<char,int> tmp = make_pair(it->nuc,it->bin);
            nucbinMap.insert(make_pair(tmp,0)).first->second++;
        }
    }
    typedef std::map<pair<char,int>,int>::const_iterator count_iter;
    for(count_iter it = nucbinMap.begin(); it != nucbinMap.end(); ++it) {
        countVec.push_back(it->second);
        nucVec.push_back(it->first.first);
        binVec.push_back(it->first.second);
    }
}

void ResultMgr::extractFromPosCache() {
    std::set<char> nucs = posCache->passingNucs(min_nuc_depth);
    int linearLengthBefore = countVec.size();
    // ABC
    if(!hasNucleotides && !hasStrands && !hasBins) // distinguish nothing
        doExtractFromPosCache<false,false,false>(nucs);
    else if(hasNucleotides && hasStrands && hasBins) // ABC
        doExtractFromPosCache<true,true,true>(nucs);
    else if(hasNucleotides && !hasStrands && !hasBins) // A
        doExtractFromPosCache<true,false,false>(nucs);
    else if(!hasNucleotides && hasStrands && !hasBins) // B
        doExtractFromPosCache<false,true,false>(nucs);
    else if(!hasNucleotides && !hasStrands && hasBins) // C
        doExtractFromPosCache<false,false,true>(nucs);
    else if(hasNucleotides && hasStrands && !hasBins) // AB
        doExtractFromPosCache<true,true,false>(nucs);
    else if(!hasNucleotides && hasStrands && hasBins) // BC
        doExtractFromPosCache<false,true,true>(nucs);
    else // AC
        doExtractFromPosCache<true,false,true>(nucs);
    int linearLengthDiff = countVec.size() - linearLengthBefore;
    // insert pos and tid if applicable
    if(linearLengthDiff > 0) {
        GenomicPosition& curGP = posCache->genomicPosition;
        posVec.insert(posVec.end(), linearLengthDiff, curGP.pos);
        if(!isRanged) {
            seqnmsVec.insert(seqnmsVec.end(), linearLengthDiff, curGP.tid + 1);
        }
    }
}

bool ResultMgr::posCachePassesFilters(const PosCache& thePosCache) {
    int totalNucFreq = thePosCache.totalNucFreq();
    int primaryNucFreq = thePosCache.primaryNucFreq();
    // must always check min_minor_allele_depth
    bool passes = totalNucFreq - primaryNucFreq >= min_minor_allele_depth;
    return passes;
}

void ResultMgr::signalGenomicPosEnd() {
    //Rprintf("start of signalPosEnd\n");
    //posCache->print();
    if(!isBuffered) {
        if(posCachePassesFilters(*posCache))
            extractFromPosCache();
        delete posCache;
        posCache = (PosCache *) NULL;
    } else {
        (*posCacheCollptrptr)->storePosCache(posCache);
        posCache = (PosCache *) NULL;
    }

    //Rprintf("end of signalPosEnd\n\n****************\n\n");
}

// extract from PosCaches corresponding to completed positions, i.e.,
// PosCache->genomicPosition < minLeftmostGenPOS
void ResultMgr::signalYieldStart() {
    //Rprintf("signalYieldStart\n");
    if(isBuffered && *posCacheCollptrptr != (PosCacheColl *) NULL) {
        for(;;) {
            posCache =
                (*posCacheCollptrptr)->destructiveNextLT(lastLeftmostGenPOS);
            if(posCache == (PosCache *) NULL)
                break;
            if(posCachePassesFilters(*posCache))
                extractFromPosCache();
            delete posCache;
            posCache = (PosCache *) NULL;
        }
    }
}

int ResultMgr::numYieldablePosCaches() const {
    //Rprintf("begin numYieldablePosCaches\n");
    if(isBuffered && (*posCacheCollptrptr != (PosCacheColl *) NULL)) {
        //Rprintf("printing posCacheColl contents\n");
        //(*posCacheCollptrptr)->printGenPositions();
        int num = (*posCacheCollptrptr)->numPosCachesLT(lastLeftmostGenPOS);
        return num;
    }
    return -1; // FIX ME: needed for end-of-input?
}

int ResultMgr::size() const {
    // number of tuples in linear representation
    return posVec.size();
}

void ResultMgr::printVecs() const {
    Rprintf("vec contents:\n");
    for(unsigned int i = 0; i != posVec.size(); ++i) {
        Rprintf("pos %d ", posVec.at(i));
        if(hasNucleotides)
            Rprintf(" nuc %c ", nucVec.at(i));
        if(hasStrands)
            Rprintf(" str %c ", strandVec.at(i));
        if(hasBins)
            Rprintf(" bin %u ", binVec.at(i));
        Rprintf(" count %d\n", countVec.at(i));
    }    
}

void ResultMgr::signalYieldEnd() {
#ifdef PILEUP_DEBUG
    printVecs();
#endif // PILEUP_DEBUG
    seqnmsVec.clear();
    posVec.clear();
    countVec.clear();
    binVec.clear();
    strandVec.clear();
    nucVec.clear();
}

void ResultMgr::signalEOI() {
    //Rprintf("got EOI message!\n");
    if(isBuffered && (*posCacheCollptrptr != (PosCacheColl *) NULL)) {
        //(*posCacheCollptrptr)->printGenPositions();
        int i = 0;
        for(;;) {
            //printf("i %d\n", i);
            posCache = (*posCacheCollptrptr)->destructiveNext();
            if(posCache == (PosCache *) NULL)
                break;
            if(posCachePassesFilters(*posCache))
                extractFromPosCache();
            delete posCache;
            posCache = (PosCache *) NULL;
            ++i;
        }
        //Rprintf("deallocating *posCacheCollptrptr's object\n");
        delete *posCacheCollptrptr;
        *posCacheCollptrptr = (PosCacheColl *) NULL;
    }
}

inline int_const_it ResultMgr::seqnmsBeg() const { return seqnmsVec.begin(); }
inline int_const_it ResultMgr::seqnmsEnd() const { return seqnmsVec.end(); }
inline int_const_it ResultMgr::posBeg() const { return posVec.begin(); }
inline int_const_it ResultMgr::posEnd() const { return posVec.end(); }
inline int_const_it ResultMgr::binBeg() const { return binVec.begin(); }
inline int_const_it ResultMgr::binEnd() const { return binVec.end(); }
inline int_const_it ResultMgr::countBeg() const { return countVec.begin(); }
inline int_const_it ResultMgr::countEnd() const { return countVec.end(); }
inline char_const_it ResultMgr::strandBeg() const { return strandVec.begin(); }
inline char_const_it ResultMgr::strandEnd() const { return strandVec.end(); }
inline char_const_it ResultMgr::nucBeg() const { return nucVec.begin(); }
inline char_const_it ResultMgr::nucEnd() const { return nucVec.end(); }
