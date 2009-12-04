#include <Rdefines.h>
#include "Alignments_class.h"

#include <limits.h> /* for CHAR_BIT */
#include <stdlib.h> /* for div() */

#define BIT7_MASK (1 << (CHAR_BIT-1))

/* Compact encoding of a logical vector */
SEXP logical_as_compact_raw_vector(SEXP x)
{
	SEXP ans;
	Rbyte *ans_elt;
	int x_length, ans_length, i, j, x_elt;
	div_t q;

	x_length = LENGTH(x);
	q = div(x_length, CHAR_BIT);
	ans_length = q.quot;
	if (q.rem != 0)
		ans_length++;
	PROTECT(ans = NEW_RAW(ans_length));
	for (i = j = 0, ans_elt = RAW(ans); i < x_length; i++, j++) {
		if (j >= CHAR_BIT) {
			j = 0;
			ans_elt++;
		}
		*ans_elt <<= 1;
		x_elt = LOGICAL(x)[i];
		if (x_elt == NA_INTEGER) {
			UNPROTECT(1);
			error("'x' contains NAs");
		}
		if (x_elt)
			*ans_elt |= 1;
	}
	if (q.rem != 0)
		*ans_elt <<= CHAR_BIT - q.rem;
	UNPROTECT(1);
	return ans;
}

/* Decoding of a compactly encoded logical vector */
SEXP compact_raw_vector_as_logical(SEXP x, SEXP length_out)
{
	SEXP ans;
	Rbyte x_elt;
	int ans_length, x_length, i, j, k;

	ans_length = INTEGER(length_out)[0];
	x_length = LENGTH(x);
	if (ans_length > x_length * CHAR_BIT)
		error("'length_out' is > 'length(x)' * %d", CHAR_BIT);
	PROTECT(ans = NEW_LOGICAL(ans_length));
	for (i = j = 0, x_elt = RAW(x)[k = 0]; i < ans_length; i++, j++) {
		if (j >= CHAR_BIT) {
			j = 0;
			x_elt = RAW(x)[++k];
		}
		LOGICAL(ans)[i] = (x_elt & BIT7_MASK) != 0;
		x_elt <<= 1;
	}
	UNPROTECT(1);
	return ans;
}

/* Subsetting of a compactly encoded logical vector */
SEXP subset_compact_raw_vector(SEXP x, SEXP subscript)
{
	SEXP ans;
	Rbyte *ans_elt;
	int x_length, subscript_length, ans_length, i, j, sub_i;
	div_t q, q2;

	x_length = LENGTH(x);
	subscript_length = LENGTH(subscript);
	q = div(subscript_length, CHAR_BIT);
	ans_length = q.quot;
	if (q.rem != 0)
		ans_length++;
	PROTECT(ans = NEW_RAW(ans_length));
	for (i = j = 0, ans_elt = RAW(ans); i < subscript_length; i++, j++) {
		if (j >= CHAR_BIT) {
			j = 0;
			ans_elt++;
		}
		*ans_elt <<= 1;
		sub_i = INTEGER(subscript)[i];
		if (sub_i == NA_INTEGER) {
			UNPROTECT(1);
			error("subscript contains NAs");
		}
		sub_i--;
		q2 = div(sub_i, CHAR_BIT);
		if (sub_i < 0 || q2.quot >= x_length) {
			UNPROTECT(1);
			error("subscript out of bounds");
		}
		if (RAW(x)[q2.quot] & (BIT7_MASK >> q2.rem))
			*ans_elt |= 1;
	}
	if (q.rem != 0)
		*ans_elt <<= CHAR_BIT - q.rem;
	UNPROTECT(1);
	return ans;
}

