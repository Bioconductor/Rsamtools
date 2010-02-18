#ifndef _ALIGNMENTS_CLASS_H_
#define _ALIGNMENTS_CLASS_H_

#include <Rinternals.h>

SEXP logical_as_compact_raw_vector(SEXP x);
SEXP compact_raw_vector_as_logical(SEXP x, SEXP length_out);
SEXP subset_compact_raw_vector(SEXP x, SEXP subscript);
SEXP compact_raw_vector_bit_count(SEXP x);
SEXP compact_raw_vector_last_bit(SEXP x);
SEXP compact_raw_vector_set_op(SEXP query, SEXP ref, SEXP align);

#endif /* _ALIGNMENTS_CLASS_H_ */
