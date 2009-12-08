#ifndef _CIGAR_H_
#define _CIGAR_H_

#include <Rinternals.h>

SEXP cigar_run_count(SEXP cigar);
SEXP cigar_table(SEXP cigar);

SEXP valid_cigar(SEXP cigar, SEXP ans_type);
SEXP split_cigar(SEXP cigar);
SEXP cigar_to_qwidth(SEXP cigar, SEXP before_hard_clipping);
SEXP cigar_to_IRanges(SEXP cigar, SEXP drop_D_ranges, SEXP merge_ranges);
SEXP cigar_to_GappedRanges(SEXP cigar, SEXP pos, SEXP flag, SEXP drop_D_ranges);
SEXP cigar_to_list_of_IRanges(SEXP cigar, SEXP rname, SEXP pos,
		SEXP flag, SEXP drop_D_ranges, SEXP merge_ranges);

#endif /* _CIGAR_H_ */
