#ifndef SCAN_BAM_DATA_H
#define SCAN_BAM_DATA_H

#include "samtools/khash.h"
#include "Rdefines.h"
#include "bam_data.h"

KHASH_SET_INIT_STR(str)

typedef struct {
    int *flag, *rname, *strand, *pos, *qwidth, *mapq, *mrnm, *mpos, *isize,
        *partition, *mates;
    char **qname;
    const char **cigar, **seq, **qual;
    khash_t(str) * cigarhash;
    int icnt, ncnt, yicnt, yncnt;
    SEXP result;
} _SCAN_BAM_DATA, *SCAN_BAM_DATA;


SCAN_BAM_DATA _Calloc_SCAN_BAM_DATA(SEXP result);
void _Free_SCAN_BAM_DATA(SCAN_BAM_DATA sbd);
int _grow_SCAN_BAM_DATA(BAM_DATA bd, int len);
void _set_mate_SCAN_BAM_DATA(int partition, int mates, void *data);
void _finish1range_SCAN_BAM_DATA(SCAN_BAM_DATA sbd, bam_header_t *header,
				 int irange);

SEXP _scan_bam_result_init(SEXP template_list, SEXP names, SEXP space);
SEXP _get_or_grow_SCAN_BAM_DATA(BAM_DATA bd, int len);

#endif
