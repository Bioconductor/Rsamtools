#ifndef PBUFFER_WRAPPER_H
#define PBUFFER_WRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif

// expose PosCacheColl dtor to C code
void pileup_pbuffer_destroy(void *pbuffer);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* PBUFFER_WRAPPER_H */
