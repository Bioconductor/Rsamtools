#ifndef PBUFFER_WRAPPER_H
#define PBUFFER_WRAPPER_H

#ifdef __cplusplus

// expose PosCacheColl dtor to C code
extern "C" void pileup_pbuffer_destroy(void *pbuffer);

#endif /* __cplusplus */
#endif /* PBUFFER_WRAPPER_H */
