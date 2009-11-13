#ifndef _CIGAR_H_
#define _CIGAR_H_

#include <Rinternals.h>

SEXP cigar_run_count(SEXP cigar);
SEXP cigar_table(SEXP cigar);
SEXP cigar_to_read_width(SEXP cigar);
SEXP cigar_to_IRanges(SEXP cigar, SEXP drop_D_ranges);
SEXP cigar_to_list_of_IRanges(SEXP cigar, SEXP rname, SEXP strand, SEXP pos,
		SEXP flag, SEXP drop_D_ranges);

#endif /* _CIGAR_H_ */
