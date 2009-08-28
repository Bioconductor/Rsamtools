#include <Rdefines.h>
#include "utilities.h"

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
    SEXP nmspc = PROTECT(_get_namespace("ShortRead"));
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
_reverse(char *buf, int len)
{
    int i;
    char tmp;
    for (i = 0; i < floor(len / 2); ++i) {
        tmp = buf[len-i-1];
        buf[len-i-1] = buf[i];
        buf[i] = tmp;
    }
}

void
_reverseComplement(char *buf, int len)
{
    static const int MAX_MAP = 256;
    static char map[256];
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
