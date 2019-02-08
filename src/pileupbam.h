#ifndef MPILEUPBAM_H
#define MPILEUPBAM_H

#include <Rdefines.h>

SEXP apply_pileups(SEXP files, SEXP names, SEXP regions, SEXP param,
                   SEXP callback);

#endif
