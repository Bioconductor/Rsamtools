#ifndef _VCFFILE_H_
#define _VCFFILE_H_

#include <Rdefines.h>

SEXP scan_vcf(SEXP ext, SEXP space, SEXP yieldSize, SEXP sample, 
              SEXP imap, SEXP gmap);
SEXP scan_vcf_connection(SEXP txt, SEXP sample, SEXP imap, SEXP gmap);

#endif                          /* _VCFFILE_H_ */
