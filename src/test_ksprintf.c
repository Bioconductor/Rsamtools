#include "test_ksprintf.h"
#include <htslib/kstring.h>

SEXP test_ksprintf(SEXP s)
{
	const char *key, *value;
	kstring_t str = {0, 0, NULL};

	key = "bcftools_viewCommand";
	value = CHAR(STRING_ELT(s, 0));
	ksprintf(&str,"##%s=%s\n", key, value);
	printf("str.s = %s\n", str.s);
	return R_NilValue;
}

