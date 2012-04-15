#ifndef _TABIXFILE_H
#define _TABIXFILE_H

#include <Rdefines.h>
#include "tabix/tabix.h"

typedef struct {
    tabix_t *tabix;
    ti_iter_t iter;
} _TABIX_FILE;

#define TABIXFILE(b) ((_TABIX_FILE *) R_ExternalPtrAddr(b))

typedef SEXP SCAN_FUN(tabix_t *tabix, ti_iter_t iter, const int size,
                      SEXP state);

SCAN_FUN tabix_as_character;

SEXP tabixfile_init();
SEXP tabixfile_open(SEXP filename, SEXP indexname);
SEXP tabixfile_close(SEXP ext);
SEXP tabixfile_isopen(SEXP ext);

SEXP index_tabix(SEXP filename, SEXP format,
                 SEXP seq, SEXP begin, SEXP end,
                 SEXP skip, SEXP comment, SEXP zeroBased);
SEXP header_tabix(SEXP ext);
SEXP scan_tabix(SEXP ext, SEXP space, SEXP yield, SEXP fun, SEXP state);

#endif
