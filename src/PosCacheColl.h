#ifndef POS_CACHE_COLL_H
#define POS_CACHE_COLL_H

#include <cstdlib>
#include <algorithm>
#include "PosCache.h"
#include <R.h>
// because set of pointers, need to compare *dereferenced* values,
// otherwise would compare pointer addresses!
struct PosCachePtrLess {
    bool operator()(const PosCache* lhs, const PosCache* rhs) const {
        return lhs->genomicPosition < rhs->genomicPosition;
    }
};

class PosCacheColl {
private:
    std::set<PosCache*, PosCachePtrLess> posCaches;
public:
    typedef std::set<PosCache*>::const_iterator const_iter;
    typedef std::set<PosCache*>::iterator iter;
    PosCacheColl() : posCaches() { }
    ~PosCacheColl() {
        while(!posCaches.empty()) {
            iter it = posCaches.begin();
            PosCache* tmp = *it;
            posCaches.erase(it);
            delete tmp;
        }
    }
#ifdef PILEUP_DEBUG
    void printGenPositions() const {
        for(const_iter it = posCaches.begin(); it != posCaches.end(); ++it) {
            printf("tid %d pos %d\n", (*it)->genomicPosition.tid,
                   (*it)->genomicPosition.pos);
        }
    }
#endif // PILEUP_DEBUG
    void storePosCache(PosCache *cachePtr) {
        //printf("storePosCache size %d\n", posCaches.size());
        if(posCaches.find(cachePtr) != posCaches.end()) {
            Rf_error("internal: posCache already in PosCacheColl");
        }
        posCaches.insert(cachePtr);
        cachePtr = (PosCache *) NULL;
    }
    // precondition: val points to an already-allocated PosCache;
    // clients will use this to check if a PosCache with matching
    // GenomicPosition already exists. If it already exists, PosCache
    // pointed to by incoming pointer will be deallocated and incoming
    // pointer will be set to existing PosCache
    // postcondition: if val found in set, remove set element
    // containing pointer, return pointer to PosCache; if not found in
    // set, return NULL
    PosCache* fetchPosCache(PosCache* val) {
        iter it = posCaches.find(val);
        if(it == posCaches.end()) {
            return val; // not found in set, leave ptr unchanged
        }
        PosCache* tmp = *it; // tmp is pointer to heap memory of PosCache
        posCaches.erase(it); // remove set<PosCache*> element
        return tmp;
    }
    PosCache* destructiveNext() {
        //printf("destructiveNext size %d\n", posCaches.size());
        if(posCaches.empty()) {
            return NULL;
        } else {
            iter it = posCaches.begin();
            PosCache* tmp = *it;
            posCaches.erase(it);
            return tmp;
        }
    }
    PosCache* destructiveNextLT(const GenomicPosition& gp) {
        if(posCaches.empty()) {
            return NULL;
        } else {
            iter it = posCaches.begin();
            if((*it)->genomicPosition < gp) {
                PosCache* tmp = *it;
                posCaches.erase(it);
                return tmp;
            } else {
                return NULL;
            }
        }
    }
    // FIX ME: use STL algorithms!!
    int numPosCachesLT(const GenomicPosition& gp) const {
        //printf("posCaches.size() %lu\n", posCaches.size());
        //printf("gp gp.tid %d gp.pos %d\n", gp.tid, gp.pos);
        if(posCaches.empty())
            return 0;
        int count = 0;
        for(const_iter it = posCaches.begin(); it != posCaches.end(); ++it) {
            if((*it)->genomicPosition < gp)
                ++count;
            else
                break;
        }
        //printf("COUNT: %d\n", count);
        return count;
    }
};

// Destructive release of ownership of PosCache in PosCacheColl that
// compares equal to posCachePtr's object
// preconditions:
// - posCachePtr points to valid initialized PosCache object
// postconditions:
// - if PosCache found in PosCacheColl:
//   a) the incoming posCachePtr's object is deallocated
//   b) posCachePtr is reassigned to found PosCache
//   c) posCacheColl no longer contains element pointing to returned
//      PosCache
// - if PosCache NOT found in PosCacheColl
//   a) posCachePtr unchanged

void getPosCacheFromColl(PosCacheColl& pcc, PosCache*& posCachePtr);

#endif /* POS_CACHE_COLL_H */
