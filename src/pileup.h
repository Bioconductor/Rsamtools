#ifndef PILEUP_H
#define PILEUP_H

#ifdef __cplusplus
extern "C" {
#endif
    SEXP c_Pileup(SEXP ext, SEXP space, SEXP keepFlags,
                  SEXP isSimpleCigar, SEXP reverseComplement,
                  SEXP yieldSize, SEXP obeyQname, SEXP asMates,
                  SEXP schema, SEXP pileupParams, SEXP seqnamesLevels);
#ifdef __cplusplus
}
#endif

#endif        
