#ifndef BAMBUFFER_H
#define BAMBUFFER_H

#include <Rinternals.h>
#include "samtools/bam.h"

typedef struct {
    bam1_t **buffer;
    int *mates, *partition;
    int i, n, as_mates, mate_flag, partition_id;
} _BAM_BUFFER, *BAM_BUFFER;

#define BAMBUFFER(b) ((BAM_BUFFER) R_ExternalPtrAddr(b))
SEXP bambuffer_init();
SEXP bambuffer(int yieldSize, int as_mates);
SEXP bambuffer_length(SEXP bufext);
SEXP bambuffer_parse(SEXP bamext, SEXP space, SEXP keepFlags, SEXP isSimpleCigar,
                     SEXP bufext, SEXP reverseComplement, SEXP template_list);
SEXP bambuffer_write(SEXP bufext, SEXP bamext, SEXP filter);

BAM_BUFFER bambuffer_new(int n, int as_mates);
void bambuffer_push(BAM_BUFFER buf, const bam1_t *bam);
void bambuffer_reset(BAM_BUFFER buf);
void bambuffer_free(BAM_BUFFER buf);

#endif
