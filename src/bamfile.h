#ifndef BAMFILE_H
#define BAMFILE_H

#include <Rdefines.h>
#include "samtools/sam.h"

typedef struct {
    samfile_t *file;
    bam_index_t *index;
    uint64_t pos0;
} _BAM_FILE;

#define BAMFILE(b) ((_BAM_FILE *) R_ExternalPtrAddr(b))

SEXP bamfile_init();
SEXP bamfile_open(SEXP filename, SEXP indexname, SEXP mode);
SEXP bamfile_close(SEXP ext);
SEXP bamfile_isopen(SEXP ext);

SEXP read_bamfile_header(SEXP ext);
SEXP scan_bamfile(SEXP ext, SEXP space, SEXP keepFlags,
                  SEXP simpleCigar, SEXP reverseComplement,
                  SEXP yieldSize, SEXP tmpl);
SEXP count_bamfile(SEXP ext, SEXP space, SEXP keepFlags, SEXP isSimpleCigar);
SEXP filter_bamfile(SEXP ext, SEXP space, SEXP keepFlags,
                    SEXP isSimpleCigar, SEXP fout_name, SEXP fout_mode);

void _check_isbamfile(SEXP ext, const char *lbl);
samfile_t *_bam_tryopen(const char *filename, const char *mode, void *aux);

#endif
