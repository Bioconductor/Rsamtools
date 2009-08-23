#ifndef _IO_SAM_H_
#define _IO_SAM_H_

#include <Rinternals.h>

/* io_sam.c */
SEXP read_bam_header(SEXP fname, SEXP mode, SEXP verbose);
SEXP scan_bam_open(SEXP fname, SEXP mode);
SEXP scan_bam_template();
SEXP scan_bam(SEXP bfile, SEXP template_list, SEXP space, 
			  SEXP keepFlags, SEXP isSimpleCigar);
SEXP count_bam(SEXP bfile, SEXP space, SEXP keepFlags, 
			   SEXP isSimpleCigar);
void scan_bam_cleanup();		/* error handling only */

#endif /* _IO_SAM_H_ */
