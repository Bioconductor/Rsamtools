#include <Rdefines.h>
#include "Alignments_class.h"

#include <limits.h> /* for CHAR_BIT */
#include <stdlib.h> /* for div() */

SEXP logical_as_compact_raw_vector(SEXP x)
{
	SEXP ans;
	Rbyte *ans_elt;
	int x_length, ans_length, i, j, k;
	div_t q;

	x_length = LENGTH(x);
	q = div(x_length, CHAR_BIT);
	ans_length = q.quot;
	if (q.rem != 0)
		ans_length++;
	PROTECT(ans = NEW_RAW(ans_length));
	for (i = j = 0, ans_elt = RAW(ans); i < ans_length; i++, ans_elt++) {
		for (k = 0; k < CHAR_BIT; k++, j++) {
			*ans_elt <<= 1;
			if (LOGICAL(x)[j])
				(*ans_elt)++;
		}
	}
	UNPROTECT(1);
	return ans;
}

