#ifndef FAFILE_H
#define FAFILE_H

#include <Rdefines.h>
#include "samtools/faidx.h"

typedef struct {
    faidx_t *index;
} _FA_FILE;

#define FAFILE(f) ((_FA_FILE *) R_ExternalPtrAddr(f))

SEXP fafile_init();
SEXP fafile_open(SEXP filename);
SEXP fafile_close(SEXP ext);
/* SEXP fafile_reopen(SEXP ext, SEXP filename, SEXP indexname); */
SEXP fafile_isopen(SEXP ext);

SEXP index_fa(SEXP filename);
SEXP n_fa(SEXP ext);
SEXP scan_fa(SEXP ext, SEXP seq, SEXP start, SEXP end, SEXP lkup);

#endif
