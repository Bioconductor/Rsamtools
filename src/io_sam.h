#ifndef _IO_SAM_H_
#define _IO_SAM_H_

#include <Rinternals.h>

/* io_sam.c */
SEXP read_bam_header(SEXP fname, SEXP mode);
SEXP scan_bam_template(SEXP tags);
SEXP index_bam(SEXP fname);
SEXP sort_bam(SEXP fname, SEXP destinationPrefix, SEXP isByQname, 
	      SEXP maxMemory);
void scan_bam_cleanup();                /* error handling only */

void _bam_check_template_list(SEXP template_list);
SEXP _read_bam_header(SEXP ext);
SEXP _scan_bam(SEXP bfile, SEXP space, SEXP keepFlags,
	       SEXP isSimpleCigar, 
	       SEXP filename, SEXP indexname, SEXP filemode,
	       SEXP reverseComplement, SEXP template_list);
SEXP _count_bam(SEXP bfile, SEXP space, SEXP keepFlags, 
		SEXP isSimpleCigar);
SEXP _filter_bam(SEXP bfile, SEXP space, SEXP keepFlags, 
		 SEXP isSimpleCigar, SEXP fout_name, SEXP fout_mode);


#endif /* _IO_SAM_H_ */
