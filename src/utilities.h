#ifndef _UTILITIES_H_
#define _UTILITIES_H_

SEXP _get_namespace(const char *pkg);
SEXP _get_strand_levels();
SEXP _get_encoding_lookup(const char *from, const char *to);
SEXP _get_lkup(const char *baseclass);
void _as_factor_SEXP(SEXP vec, SEXP lvls);
void _as_factor(SEXP vec, const char **lvls, const int n_lvls);
SEXP _as_XStringSet(const char **key, int len, const char *baseclass);
SEXP _as_PhredQuality(const char **key, int len);
void _reverse(char *buf, int len);
void _reverseComplement(char *buf, int len);
char *_rtrim(char *);

/* common checks */

void _checkext(SEXP ext, SEXP tag, const char *lbl);
void _checknames(SEXP filename, SEXP indexname, SEXP filemode);
void _checkparams(SEXP space, SEXP keepFlags, SEXP isSimpleCigar);

/* pairing */

SEXP p_pairing(SEXP x_qname, SEXP x_flag, SEXP x_rname,
               SEXP x_pos, SEXP x_rnext, SEXP x_pnext,
               SEXP y_qname, SEXP y_flag, SEXP y_rname,
               SEXP y_pos, SEXP y_rnext, SEXP y_pnext);

SEXP find_mate_within_groups(SEXP group_sizes,
                             SEXP x_flag, SEXP x_rname,
                             SEXP x_pos, SEXP x_rnext, SEXP x_pnext);

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

#endif                          /* _UTILITIES_H_ */
