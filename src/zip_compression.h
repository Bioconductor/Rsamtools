#ifndef ZIP_COMPRESSION_H
#define ZIP_COMPRESSION_H

#include <Rdefines.h>

SEXP bgzip(SEXP from, SEXP dest);
SEXP razip(SEXP from, SEXP dest);

#endif
