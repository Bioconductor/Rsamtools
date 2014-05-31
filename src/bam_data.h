#ifndef BAM_DATA_H
#define BAM_DATA_H

#include "Rdefines.h"
#include "bamfile.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int BLOCKSIZE;              /* size to grow vectors */
    char *CIGAR_BUF;            /* string representation of CIGAR */
    uint32_t CIGAR_BUF_SZ;

    int parse_status;
    BAM_FILE bfile;
    int irec, iparsed, irange, nrange;
    uint32_t keep_flag[2], cigar_flag;
    int reverseComplement, yieldSize, 
        obeyQname, asMates;
    char qnamePrefixEnd, qnameSuffixStart;

    void *extra;
} _BAM_DATA, *BAM_DATA;

enum {
    QNAME_IDX = 0, FLAG_IDX, RNAME_IDX, STRAND_IDX, POS_IDX, QWIDTH_IDX,
    MAPQ_IDX, CIGAR_IDX, MRNM_IDX, MPOS_IDX, ISIZE_IDX, SEQ_IDX,
    QUAL_IDX, TAG_IDX, PARTITION_IDX, MATES_IDX
};

BAM_DATA _init_BAM_DATA(SEXP ext, SEXP space, SEXP flag, SEXP isSimpleCigar,
			int reverseComplement, int yieldSize, int obeyQname,
			int asMates, char qnamePrefixEnd, 
                        char qnameSuffixStart, void *extra);
void _Free_BAM_DATA(BAM_DATA bd);
BAM_FILE _bam_file_BAM_DATA(BAM_DATA bd);
int _count1_BAM_DATA(const bam1_t *bam, BAM_DATA bd);
int _filter_and_parse1_BAM_DATA(const bam1_t *bam, BAM_DATA bd);
int _filter1_BAM_DATA(const bam1_t *bam, BAM_DATA bd);
int _parse1_BAM_DATA(const bam1_t *bam, BAM_DATA bd);
void _finish1range_BAM_DATA(BAM_DATA  bd);

#ifdef __cplusplus
}
#endif

#endif
