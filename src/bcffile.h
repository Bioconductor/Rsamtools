#ifndef _BCFFILE_H_
#define _BCFFILE_H_

#include <Rdefines.h>
#include "bcftools/bcf.h"

/* io_bcf.c */
typedef struct {
    bcf_t *file;
    bcf_idx_t *index;
} _BCF_FILE;

#define BCFFILE(b) ((_BCF_FILE *) R_ExternalPtrAddr(b))

SEXP bcffile_init();
SEXP bcffile_open(SEXP filename, SEXP indexname, SEXP mode);
SEXP bcffile_close(SEXP ext);
SEXP bcffile_isopen(SEXP ext);
SEXP bcffile_isvcf(SEXP ext);

SEXP scan_bcf_header(SEXP ext);
SEXP scan_bcf(SEXP ext, SEXP space, SEXP typemap);

SEXP as_bcf(SEXP file, SEXP dictionary, SEXP destination);
SEXP index_bcf(SEXP file);

#endif                          /* _BCFFILE_H_ */
