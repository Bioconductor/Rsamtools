#include <Rdefines.h>
#include "utilities.h"

void
_samtools_exit(int status)
{
    Rf_error("internal: samtools invoked 'exit(%d)'; see warnings()", 
	     status);
}

int
_samtools_fprintf(FILE *file, const char *fmt, ...)
{
    static const int bufsize = 2048;
    va_list argp;
    int n;

    if (stderr != file) {
	va_start(argp, fmt);
	n = vfprintf(file, fmt, argp);
	va_end(argp);
    } else {
	/* silence some messages */
	char *buf = (char *) R_alloc(bufsize, sizeof(char));
	if (0 == strncmp("[samopen] SAM header is present:", fmt, 32))
	    return 0;
	va_start(argp, fmt);
	n = vsnprintf(buf, bufsize, fmt, argp);
	va_end(argp);
	Rf_warning(buf);
    }
    return n;
}

SEXP
_get_namespace(const char* pkg)
{
    SEXP fun = PROTECT(findFun(install("getNamespace"), R_GlobalEnv));
    SEXP nmspc = PROTECT(NEW_CHARACTER(1));
    SET_STRING_ELT(nmspc, 0, mkChar(pkg));
    nmspc = eval(lang2(fun, nmspc), R_GlobalEnv);
    UNPROTECT(2);
    return nmspc;
}


SEXP
_get_strand_levels()
{
    SEXP nmspc = PROTECT(_get_namespace("Rsamtools"));
    SEXP ans = eval(findVar(install(".STRAND_LEVELS"), nmspc), nmspc);
    UNPROTECT(1);
    return ans;
}

void
_as_factor_SEXP(SEXP vec, SEXP lvls)
{
    SEXP cls = PROTECT(NEW_CHARACTER(1));
    SET_STRING_ELT(cls, 0, mkChar("factor"));
    SET_CLASS(vec, cls);
    SET_ATTR(vec, install("levels"), lvls);
    UNPROTECT(1);
}

void
_as_factor(SEXP vec, const char **lvls, const int n_lvls)
{
    SEXP levels = PROTECT(NEW_CHARACTER(n_lvls));
    int i;
    for (i = 0; i < n_lvls; ++i)
        SET_STRING_ELT(levels, i, mkChar(lvls[i]));
    _as_factor_SEXP(vec, levels);
    UNPROTECT(1);
}

void
_reverse(unsigned char *buf, int len)
{
    int i;
    unsigned char tmp;
    for (i = 0; i < floor(len / 2); ++i) {
        tmp = buf[len-i-1];
        buf[len-i-1] = buf[i];
        buf[i] = tmp;
    }
}

void
_reverseComplement(unsigned char *buf, int len)
{
    static const int MAX_MAP = 256;
    static unsigned char map[256];
    static int init = 0;
    if (init == 0) {
        init = 1;
        for (int i = 0; i < MAX_MAP; ++i)
            map[i] = (char) i;
        map['A'] = 'T'; map['C'] = 'G'; map['G'] = 'C'; map['T'] = 'A';
        map['a'] = 't'; map['c'] = 'g'; map['g'] = 'c'; map['t'] = 'a';
    }
    _reverse(buf, len);
    for (int i = 0; i < len; ++i)
        buf[i] = map[(int) buf[i]];
}

/* scan-related */

void
_scan_checkext(SEXP ext, SEXP tag, const char *lbl)
{
    if (EXTPTRSXP != TYPEOF(ext) || 
        tag != R_ExternalPtrTag(ext))
        Rf_error("incorrect instance for '%s'", lbl);
}

void
_scan_checknames(SEXP filename, SEXP indexname, SEXP filemode)
{
    if (!IS_CHARACTER(filename) || LENGTH(filename) > 1)
        Rf_error("'filename' must be character(0) or character(1)");
    if (!IS_CHARACTER(indexname) || LENGTH(indexname) > 1)
        Rf_error("'indexname' must be character(0) or character(1)");
    if (!IS_CHARACTER(filemode) || LENGTH(filemode) != 1)
        Rf_error("'filemode' must be character(1)");
}

void
_scan_checkparams(SEXP space, SEXP keepFlags, SEXP isSimpleCigar)
{
    const int MAX_CHRLEN = 1 << 29; /* See samtools/bam_index.c */
    if (R_NilValue != space) {
        if (!IS_VECTOR(space) || LENGTH(space) != 3)
            Rf_error("'space' must be list(3) or NULL");
        if (!IS_CHARACTER(VECTOR_ELT(space, 0)))
            Rf_error("internal: 'space[1]' must be character()");
        if (!IS_INTEGER(VECTOR_ELT(space, 1)))
            Rf_error("internal: 'space[2]' must be integer()");
        if (!IS_INTEGER(VECTOR_ELT(space, 2)))
            Rf_error("internal: 'space[3]' must be integer()");
        if ((LENGTH(VECTOR_ELT(space, 0)) != LENGTH(VECTOR_ELT(space, 1))) ||
            (LENGTH(VECTOR_ELT(space, 0)) != LENGTH(VECTOR_ELT(space, 2))))
            Rf_error("internal: 'space' elements must all be the same length");
        const int *end = INTEGER(VECTOR_ELT(space, 2)),
            nrange = LENGTH(VECTOR_ELT(space, 2));
        for (int irange = 0; irange < nrange; ++irange)
            if (end[irange] > MAX_CHRLEN)
                Rf_error("'end' must be <= %d", MAX_CHRLEN);
    }
    if (R_NilValue != keepFlags)
	if (!IS_INTEGER(keepFlags) || LENGTH(keepFlags) != 2)
	    Rf_error("'keepFlags' must be integer(2) or NULL");
    if (R_NilValue != isSimpleCigar)
	if (!IS_LOGICAL(isSimpleCigar)  || LENGTH(isSimpleCigar) != 1)
	    Rf_error("'isSimpleCigar' must be logical(1) or NULL");
}
