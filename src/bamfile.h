#ifndef BAMFILE_H
#define BAMFILE_H

#include <Rdefines.h>
#include <samtools-1.7-compat.h>
#include <htslib/hfile.h>
#include "bambuffer.h"
#include "bam_mate_iter.h"
#include "pbuffer_wrapper.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    htsFile *file;
    hts_idx_t *index;
    bam_hdr_t *header;
    int64_t pos0;
    int irange0;
    bam_mate_iter_t iter;
    void *pbuffer; /* for buffered pileup */
} _BAM_FILE, *BAM_FILE;

#define BAMFILE(b) ((BAM_FILE) R_ExternalPtrAddr(b))

SEXP bamfile_init();
SEXP bamfile_open(SEXP file0, SEXP file1, SEXP mode);
SEXP bamfile_close(SEXP ext);
SEXP bamfile_isopen(SEXP ext);
SEXP bamfile_isincomplete(SEXP ext);

SEXP read_bamfile_header(SEXP ext, SEXP what);
SEXP scan_bamfile(SEXP ext, SEXP regions, SEXP keepFlags,
                  SEXP simpleCigar, SEXP tagFilter,  SEXP mapqFilter,
                  SEXP reverseComplement, SEXP yieldSize,
                  SEXP tmpl, SEXP obeyQname,
                  SEXP asMates, SEXP qnamePrefix, SEXP qnameSuffix);
SEXP count_bamfile(SEXP ext, SEXP regions, SEXP keepFlags, SEXP isSimpleCigar,
                   SEXP tagFilter, SEXP mapqFilter);
SEXP prefilter_bamfile(SEXP ext, SEXP regions, SEXP keepFlags,
                       SEXP isSimpleCigar, SEXP tagFilter, SEXP mapqFilter,
                       SEXP yieldSize,
                       SEXP obeyQname, SEXP asMates, SEXP qnamePrefix,
                       SEXP qnameSuffix);
SEXP filter_bamfile(SEXP ext, SEXP regions, SEXP keepFlags,
                    SEXP isSimpleCigar, SEXP tagFilter, SEXP mapqFilter,
                    SEXP ext_out);

void _check_isbamfile(SEXP ext, const char *lbl);
htsFile *_bam_tryopen(const char *filename, const char *mode, void *aux);

#ifdef __cplusplus
}
#endif

#endif
