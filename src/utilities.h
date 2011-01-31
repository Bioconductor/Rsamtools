#ifndef _UTILITIES_H_
#define _UTILITIES_H_

/* capture samtools errors */
extern void _samtools_abort();
extern void _samtools_exit(int status);
extern int _samtools_fprintf (FILE *, const char *, ...);

SEXP _get_namespace(const char *pkg);
SEXP _get_strand_levels();
void _as_factor_SEXP(SEXP vec, SEXP lvls);
void _as_factor(SEXP vect, const char **lvls, const int n_lvls);
void _reverse(unsigned char *buf, int len);
void _reverseComplement(unsigned char *buf, int len);

/* scan-related */

void _scan_checkext(SEXP ext, SEXP tag, const char *lbl);
void _scan_checknames(SEXP filename, SEXP indexname, SEXP filemode);
void _scan_checkparams(SEXP space, SEXP keepFlags, SEXP isSimpleCigar);

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
