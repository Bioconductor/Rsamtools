#ifndef _UTILITIES_H_
#define _UTILITIES_H_

SEXP _get_namespace(const char *pkg);
SEXP _get_strand_levels();
void _as_factor_SEXP(SEXP vec, SEXP lvls);
void _as_factor(SEXP vect, const char **lvls, const int n_lvls);
void _reverse(char *buf, int len);
void _reverseComplement(char *buf, int len);

/* call-building macros */

#define NEW_CALL(S, T, NAME, ENV, N)            \
    PROTECT(S = T = allocList(N));              \
    SET_TYPEOF(T, LANGSXP);                     \
    SETCAR(T, findFun(install(NAME), ENV));     \
    T = CDR(T)
#define CSET_CDR(T, NAME, VALUE)                \
    SETCAR(T, VALUE);                           \
    SET_TAG(T, install(NAME));                  \
    T = CDR(T)
#define CEVAL_TO(S, ENV, GETS)                  \
    GETS = eval(S, ENV);                        \
    UNPROTECT(1)

#endif /* _UTILITIES_H_ */
