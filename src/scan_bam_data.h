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
    khash_t(str) *cigarhash;
    int icnt, ncnt,
        mates_flag, partition_id, /* set prior to parsing 1 bam record */
        partition_as_width;
    SEXP result;
} _SCAN_BAM_DATA, *SCAN_BAM_DATA;


SCAN_BAM_DATA _init_SCAN_BAM_DATA(SEXP result, int partition_as_width);
void _Free_SCAN_BAM_DATA(SCAN_BAM_DATA sbd);
int _grow_SCAN_BAM_DATA(BAM_DATA bd, int len);
void _finish1range_SCAN_BAM_DATA(SCAN_BAM_DATA sbd, bam_header_t *header,
				 int irange);

SEXP _scan_bam_result_init(SEXP template_list, SEXP names, SEXP space);
SEXP _get_or_grow_SCAN_BAM_DATA(BAM_DATA bd, int len);

#endif
