#include "scan_bam_data.h"
#include "utilities.h"

/* _SCAN_BAM_DATA */

SCAN_BAM_DATA _init_SCAN_BAM_DATA(SEXP result)
{
    SCAN_BAM_DATA sbd = Calloc(1, _SCAN_BAM_DATA);
    sbd->cigarhash = kh_init(str);
    sbd->result = result;
    sbd->mates_flag = NA_LOGICAL;
    sbd->partition_id = 0;
    return sbd;
}

static void _Free_strhash(khash_t(str) * h)
{
    khiter_t k;
    char *buf;
    for (k = kh_begin(h); kh_end(h) != k; ++k)
        if (kh_exist(h, k)) {
            buf = (char *) kh_key(h, k);
            Free(buf);
        }
    kh_destroy(str, h);
}

void _Free_SCAN_BAM_DATA(SCAN_BAM_DATA sbd)
{
    _Free_strhash(sbd->cigarhash);
    Free(sbd);
}

static void _grow_SCAN_BAM_DATA_tags(SEXP tags, int len)
{
    int i, j;

    for (i = 0; i < LENGTH(tags); ++i) {
        SEXP elt = VECTOR_ELT(tags, i);
        int len0 = LENGTH(elt);
        elt = Rf_lengthgets(elt, len);
        SET_VECTOR_ELT(tags, i, elt);
        switch (TYPEOF(elt)) {
        case INTSXP:
            for (j = len0; j < len; ++j)
                INTEGER(elt)[j] = NA_INTEGER;
            break;
        case REALSXP:
            for (j = len0; j < len; ++j)
                REAL(elt)[j] = NA_REAL;
            break;
        case STRSXP:
            for (j = len0; j < len; ++j)
                SET_STRING_ELT(elt, j, NA_STRING);
            break;
        case RAWSXP:
            for (j = len0; j < len; ++j)
                RAW(elt)[j] = 0x0;
        default:
            break;
        }
    }
}

int _grow_SCAN_BAM_DATA(BAM_DATA bd, int len)
{
    int i;
    SEXP r, s;
    SCAN_BAM_DATA sbd = (SCAN_BAM_DATA) bd->extra;

    r = VECTOR_ELT(sbd->result, bd->irange);

    for (i = 0; i < LENGTH(r); ++i) {
        if (R_NilValue == (s = VECTOR_ELT(r, i)))
            continue;
        switch (i) {
        case FLAG_IDX:
            sbd->flag = _Rs_Realloc(sbd->flag, len, int);
            break;
        case RNAME_IDX:
            sbd->rname = _Rs_Realloc(sbd->rname, len, int);
            break;
        case STRAND_IDX:
            sbd->strand = _Rs_Realloc(sbd->strand, len, int);
            break;
        case POS_IDX:
            sbd->pos = _Rs_Realloc(sbd->pos, len, int);
            break;
        case QWIDTH_IDX:
            sbd->qwidth = _Rs_Realloc(sbd->qwidth, len, int);
            break;
        case MAPQ_IDX:
            sbd->mapq = _Rs_Realloc(sbd->mapq, len, int);
            break;
        case MRNM_IDX:
            sbd->mrnm = _Rs_Realloc(sbd->mrnm, len, int);
            break;
        case MPOS_IDX:
            sbd->mpos = _Rs_Realloc(sbd->mpos, len, int);
            break;
        case ISIZE_IDX:
            sbd->isize = _Rs_Realloc(sbd->isize, len, int);
            break;
        case QNAME_IDX:
            sbd->qname = _Rs_Realloc(sbd->qname, len, char *);
            break;
        case CIGAR_IDX:
            sbd->cigar = _Rs_Realloc(sbd->cigar, len, const char *);
            break;
        case SEQ_IDX:
            sbd->seq = _Rs_Realloc(sbd->seq, len, const char *);
            break;
        case QUAL_IDX:
            sbd->qual = _Rs_Realloc(sbd->qual, len, const char *);
            break;
        case TAG_IDX:
            if (R_NilValue != s)
                _grow_SCAN_BAM_DATA_tags(s, len);
            break;
        case PARTITION_IDX:
            sbd->partition = _Rs_Realloc(sbd->partition, len, int);
            break;
        case MATES_IDX:
            sbd->mates = _Rs_Realloc(sbd->mates, len, int);
            break;
        default:
            Rf_error("[Rsamtools internal] unhandled _grow_SCAN_BAM_DATA");
            break;
        }
    }

    return len;
}

SEXP _get_or_grow_SCAN_BAM_DATA(BAM_DATA bd, int len)
{
    _SCAN_BAM_DATA *sbd = (_SCAN_BAM_DATA *) bd->extra;
    if (len < 0) {
        if (sbd->icnt < sbd->ncnt)
            return VECTOR_ELT(sbd->result, bd->irange);
        len = sbd->ncnt + bd->BLOCKSIZE;
    }

    sbd->ncnt = _grow_SCAN_BAM_DATA(bd, len);
    return VECTOR_ELT(sbd->result, bd->irange);
}

void _finish1range_SCAN_BAM_DATA(SCAN_BAM_DATA sbd, bam_header_t *header,
				 int irange)
{
    /* FIXME: replace mrnm '=' with rname */
    const char *mates_lvls[] = { "mated", "ambiguous", "unmated" };
    int i, j;
    SEXP r, s;
    r = VECTOR_ELT(sbd->result, irange);
    for (i = 0; i < LENGTH(r); ++i) {
        if (R_NilValue == (s = VECTOR_ELT(r, i)))
            continue;
        switch (i) {
        case FLAG_IDX:
            s = Rf_lengthgets(s, sbd->icnt);
            SET_VECTOR_ELT(r, i, s);
            memcpy(INTEGER(s), sbd->flag, sbd->icnt * sizeof(int));
            Free(sbd->flag);
            break;
        case RNAME_IDX:
            s = Rf_lengthgets(s, sbd->icnt);
            SET_VECTOR_ELT(r, i, s);
            memcpy(INTEGER(s), sbd->rname, sbd->icnt * sizeof(int));
            _as_rname(s, (const char **) header->target_name,
                      header->n_targets);
            Free(sbd->rname);
            break;
        case STRAND_IDX:
            s = Rf_lengthgets(s, sbd->icnt);
            SET_VECTOR_ELT(r, i, s);
            memcpy(INTEGER(s), sbd->strand, sbd->icnt * sizeof(int));
            _as_strand(s);
            Free(sbd->strand);
            break;
        case POS_IDX:
            s = Rf_lengthgets(s, sbd->icnt);
            SET_VECTOR_ELT(r, i, s);
            memcpy(INTEGER(s), sbd->pos, sbd->icnt * sizeof(int));
            Free(sbd->pos);
            break;
        case QWIDTH_IDX:
            s = Rf_lengthgets(s, sbd->icnt);
            SET_VECTOR_ELT(r, i, s);
            memcpy(INTEGER(s), sbd->qwidth, sbd->icnt * sizeof(int));
            Free(sbd->qwidth);
            break;
        case MAPQ_IDX:
            s = Rf_lengthgets(s, sbd->icnt);
            SET_VECTOR_ELT(r, i, s);
            memcpy(INTEGER(s), sbd->mapq, sbd->icnt * sizeof(int));
            Free(sbd->mapq);
            break;
        case MRNM_IDX:
            s = Rf_lengthgets(s, sbd->icnt);
            SET_VECTOR_ELT(r, i, s);
            memcpy(INTEGER(s), sbd->mrnm, sbd->icnt * sizeof(int));
            _as_rname(s, (const char **) header->target_name,
                      header->n_targets);
            Free(sbd->mrnm);
            break;
        case MPOS_IDX:
            s = Rf_lengthgets(s, sbd->icnt);
            SET_VECTOR_ELT(r, i, s);
            memcpy(INTEGER(s), sbd->mpos, sbd->icnt * sizeof(int));
            Free(sbd->mpos);
            break;
        case ISIZE_IDX:
            s = Rf_lengthgets(s, sbd->icnt);
            SET_VECTOR_ELT(r, i, s);
            memcpy(INTEGER(s), sbd->isize, sbd->icnt * sizeof(int));
            Free(sbd->isize);
            break;
        case QNAME_IDX:
            s = Rf_lengthgets(s, sbd->icnt);
            SET_VECTOR_ELT(r, i, s);
            for (j = 0; j < sbd->icnt; ++j) {
                SET_STRING_ELT(s, j, mkChar(sbd->qname[j]));
                Free(sbd->qname[j]);
            }
            Free(sbd->qname);
            break;
        case CIGAR_IDX:
            s = Rf_lengthgets(s, sbd->icnt);
            SET_VECTOR_ELT(r, i, s);
            for (j = 0; j < sbd->icnt; ++j)
                if (NULL == sbd->cigar[j])
                    SET_STRING_ELT(s, j, NA_STRING);
                else
                    SET_STRING_ELT(s, j, mkChar(sbd->cigar[j]));
            Free(sbd->cigar);
            break;
        case SEQ_IDX:
            s = _as_XStringSet(sbd->seq, sbd->icnt, "DNAString");
            SET_VECTOR_ELT(r, i, s);
            for (j = 0; j < sbd->icnt; ++j)
                Free(sbd->seq[j]);
            Free(sbd->seq);
            break;
        case QUAL_IDX:
            s = _as_PhredQuality(sbd->qual, sbd->icnt);
            SET_VECTOR_ELT(r, i, s);
            for (j = 0; j < sbd->icnt; ++j)
                Free(sbd->qual[j]);
            Free(sbd->qual);
            break;
        case TAG_IDX:
            _grow_SCAN_BAM_DATA_tags(s, sbd->icnt);
            break;
        case PARTITION_IDX:
            s = Rf_lengthgets(s, sbd->icnt);
            SET_VECTOR_ELT(r, i, s);
            memcpy(INTEGER(s), sbd->partition, Rf_length(s) * sizeof(int));
            Free(sbd->partition);
            break;
        case MATES_IDX:
            s = Rf_lengthgets(s, sbd->icnt);
            SET_VECTOR_ELT(r, i, s);
            memcpy(INTEGER(s), sbd->mates, sbd->icnt * sizeof(int));
            _as_factor(s, mates_lvls, 3);
            Free(sbd->mates);
            break;
        default:
            Rf_error("[Rsamtools internal] unhandled _finish1range_BAM_DATA");
            break;
        }
    }

    sbd->icnt = sbd->ncnt = 0;
    sbd->mates_flag = NA_LOGICAL;
}
