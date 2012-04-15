#ifndef _IO_SAM_H_
#define _IO_SAM_H_

#include <Rinternals.h>

SEXP read_bam_header(SEXP fname, SEXP mode);
SEXP scan_bam_template(SEXP tags);
SEXP as_bam(SEXP fname, SEXP destination);
SEXP sort_bam(SEXP fname, SEXP destinationPrefix, SEXP isByQname,
              SEXP maxMemory);
SEXP merge_bam(SEXP fnames, SEXP destination, SEXP overwrite,
               SEXP hname, SEXP regionStr, SEXP isByQname,
               SEXP addRG, SEXP compressLevel1);
SEXP index_bam(SEXP indexname);
void scan_bam_cleanup();        /* error handling only */

void _bam_check_template_list(SEXP template_list);
SEXP _read_bam_header(SEXP ext);
SEXP _scan_bam(SEXP bfile, SEXP space, SEXP keepFlags,
               SEXP isSimpleCigar, SEXP reverseComplement,
               SEXP yieldSize, SEXP template_list);
SEXP _count_bam(SEXP bfile, SEXP space, SEXP keepFlags, SEXP isSimpleCigar);
SEXP _filter_bam(SEXP bfile, SEXP space, SEXP keepFlags,
                 SEXP isSimpleCigar, SEXP fout_name, SEXP fout_mode);

#endif                          /* _IO_SAM_H_ */
