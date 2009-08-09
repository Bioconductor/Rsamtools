#ifndef _RSAMTOOLS_H_
#define _RSAMTOOLS_H_

#include <Rdefines.h>
#include <R_ext/Rdynload.h>

/* io_sam.c */
SEXP read_bam_header(SEXP fnaem, SEXP mode, SEXP verbose);
SEXP scan_bam_open(SEXP fname, SEXP mode);
SEXP scan_bam_template();
SEXP scan_bam(SEXP bfile, SEXP template_list, SEXP space, SEXP flag);

/* utilities.c */
SEXP _get_namespace(const char *pkg);
SEXP _get_strand_levels();
void _as_factor_SEXP(SEXP vec, SEXP lvls);
void _as_factor(SEXP vect, const char **lvls, const int n_lvls);
void _reverse(char *buf);
void _reverseComplement(char *buf);

#endif /* _RSAMTOOLS_H_ */
