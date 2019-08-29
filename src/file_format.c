#include "file_format.h"
#include <stdlib.h>
#include <htslib/hts.h>

SEXP file_format(SEXP ext)
{
    if (!IS_CHARACTER(ext) || Rf_length(ext) != 1)
        Rf_error("fileFormat 'file' must be character(1)");
    const char *fname = translateChar(STRING_ELT(ext, 0));

    htsFile *hts = hts_open(fname, "r");
    if (hts == NULL)
        Rf_error("'fileFormat()' failed to open file '%s'", fname);
    const htsFormat *format = hts_get_format(hts);
    const char *type = hts_format_description(format);

    SEXP result = PROTECT(Rf_ScalarString(mkChar(type)));

    if (hts_close(hts) != 0)
        Rf_error("'fileFormat()' failed to close file '%s'", fname);
    free((void *) type);        /* FIXME: could fail on Windows */

    UNPROTECT(1);
    return result;
}
