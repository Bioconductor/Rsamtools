#ifndef PILEUP_H
#define PILEUP_H

#ifdef __cplusplus
#define R_NO_REMAP
#include "io_sam.h"
#include "utilities.h"
#include "PileupBufferShim.h"
#ifdef PILEUP_DEBUG
#include "nate_utilities.h"
#endif

extern "C" {
#endif
    SEXP c_Pileup(SEXP ext, SEXP regions, SEXP keepFlags,
                  SEXP isSimpleCigar, SEXP tagFilter, SEXP mapqFilter,
                  SEXP reverseComplement,
                  SEXP yieldSize, SEXP obeyQname, SEXP asMates,
                  SEXP qnamePrefixEnd, SEXP qnameSuffixStart, 
                  SEXP pileupParams);
#ifdef __cplusplus
}
#endif

#endif        
