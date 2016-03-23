#include "PosCacheColl.h"

// because it's a free function we have to ensure it's only ever
// defined once. Including it in the .h triggers multiple definition
// linking error
void getPosCacheFromColl(PosCacheColl& pcc, PosCache*& posCachePtr) {
    PosCache* tmp = posCachePtr; // hold address pointed to by posCachePtr
    posCachePtr = pcc.fetchPosCache(posCachePtr);
    if(tmp != posCachePtr) {
        // found val in pcc, must deallocate incoming ptr's object
        delete tmp;
        tmp = (PosCache *) NULL;
    }
}
