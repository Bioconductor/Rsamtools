#ifndef _ALIGNMENTS_CLASS_H_
#define _ALIGNMENTS_CLASS_H_

#include <Rinternals.h>

SEXP logical_as_compact_raw_vector(SEXP x);
SEXP compact_raw_vector_as_logical(SEXP x, SEXP length_out);

#endif /* _ALIGNMENTS_CLASS_H_ */
