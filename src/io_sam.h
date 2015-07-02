#ifndef _IO_SAM_H_
#define _IO_SAM_H_

#include <Rinternals.h>
#include "bam_data.h"

#ifdef __cplusplus
extern "C" {
#endif

SEXP scan_bam_template(SEXP rname, SEXP tags);
SEXP sort_bam(SEXP fname, SEXP destinationPrefix, SEXP isByQname,
              SEXP maxMemory);
SEXP merge_bam(SEXP fnames, SEXP destination, SEXP overwrite,
               SEXP hname, SEXP regionStr, SEXP isByQname,
               SEXP addRG, SEXP compressLevel1);
SEXP index_bam(SEXP indexname);
void scan_bam_cleanup();        /* error handling only */

void _bam_check_template_list(SEXP template_list);
SEXP _read_bam_header(SEXP ext, SEXP what);
SEXP _scan_bam(SEXP bfile, SEXP space, SEXP keepFlags,
               SEXP isSimpleCigar, SEXP tagFilter, SEXP mapqFilter,
               SEXP reverseComplement, SEXP yieldSize,
               SEXP template_list, SEXP obeyQname, SEXP asMates,
               SEXP qnamePrefixEnd, SEXP qnameSuffixStart);
SEXP _count_bam(SEXP bfile, SEXP space, SEXP keepFlags, SEXP isSimpleCigar,
                SEXP tagFilter, SEXP mapqFilter);
SEXP _prefilter_bam(SEXP bfile, SEXP space, SEXP keepFlags,
		    SEXP isSimpleCigar, SEXP tagFilter, SEXP mapqFilter,
                    SEXP yieldSize, SEXP obeyQname, SEXP asMates,
                    SEXP qnamePrefixEnd, SEXP qnameSuffixStart);
SEXP _filter_bam(SEXP bfile, SEXP space, SEXP keepFlags,
                 SEXP isSimpleCigar, SEXP tagFilter, SEXP mapqFilter,
                 SEXP fout_name, SEXP fout_mode);

typedef void (_FINISH1_FUNC) (BAM_DATA);
int _do_scan_bam(BAM_DATA bd, SEXP space, bam_fetch_f parse1,
                 bam_fetch_mate_f parse1_mate, _FINISH1_FUNC finish1);


#ifdef __cplusplus
}
#endif

#endif                          /* _IO_SAM_H_ */
