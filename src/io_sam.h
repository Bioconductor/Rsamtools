#ifndef _IO_SAM_H_
#define _IO_SAM_H_

#include <Rinternals.h>

/* io_sam.c */
SEXP read_bam_header(SEXP fname, SEXP mode);
SEXP scan_bam_template(SEXP tags);
SEXP scan_bam(SEXP fname, SEXP index, SEXP mode,
			  SEXP space, SEXP keepFlags, SEXP isSimpleCigar,
			  SEXP reverseComplement, SEXP template_list);
SEXP filter_bam(SEXP fname, SEXP index, SEXP mode,
				SEXP space, SEXP keepFlags, SEXP isSimpleCigar,
				SEXP fout_name, SEXP fout_mode);
SEXP count_bam(SEXP fname, SEXP index, SEXP mode, 
			   SEXP space, SEXP keepFlags, SEXP isSimpleCigar);
SEXP index_bam(SEXP fname);
SEXP sort_bam(SEXP fname, SEXP destinationPrefix, SEXP isByQname, SEXP maxMemory);
void scan_bam_cleanup();		/* error handling only */

#endif /* _IO_SAM_H_ */
