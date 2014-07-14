#include "pbuffer_wrapper.h"
#include "PosCacheColl.h"

void pileup_pbuffer_destroy(void *pbuffer) {
    if(pbuffer != NULL) {
        delete (static_cast<PosCacheColl*>(pbuffer));
    }
}
