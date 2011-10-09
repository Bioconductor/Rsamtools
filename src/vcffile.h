#ifndef _VCFFILE_H_
#define _VCFFILE_H_

#include <Rdefines.h>

SEXP scan_vcf(SEXP ext, SEXP space, SEXP yieldSize, SEXP sample, SEXP map);
SEXP scan_vcf_connection(SEXP txt, SEXP sample, SEXP map);

#endif                          /* _VCFFILE_H_ */
