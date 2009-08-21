#ifndef _RSAMTOOLS_H_
#define _RSAMTOOLS_H_

#include <Rdefines.h>
#include <R_ext/Rdynload.h>

/* io_sam.c */
SEXP read_bam_header(SEXP fname, SEXP mode, SEXP verbose);
SEXP scan_bam_open(SEXP fname, SEXP mode);
SEXP scan_bam_template();
SEXP scan_bam(SEXP bfile, SEXP template_list, SEXP space, 
			  SEXP keepFlags, SEXP isSimpleCigar);

void scan_bam_cleanup();		/* error handling only */

#endif /* _RSAMTOOLS_H_ */
