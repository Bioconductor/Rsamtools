#ifndef MPILEUPBAM_H
#define MPILEUPBAM_H

#include <Rdefines.h>

SEXP pileup_bam(SEXP files, SEXP space, SEXP param, 
		 SEXP callback);

#endif
