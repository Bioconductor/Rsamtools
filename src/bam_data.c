#include <htslib/khash.h>
#include <sam.h>
#include "bamfile.h"
#include "bam_data.h"
#include "scan_bam_data.h"
#include "utilities.h"

#define BAM_PARSE_STATUS_OK 0

static const int BAM_INIT_SIZE = 1048576;

enum { CIGAR_SIMPLE = 1 };

/* _BAM_DATA */

static BAM_DATA _Calloc_BAM_DATA(int blocksize, int cigar_buf_sz)
{
    BAM_DATA bd = Calloc(1, _BAM_DATA);
    bd->BLOCKSIZE = blocksize;
    bd->cigar_buf_sz = cigar_buf_sz;
    bd->cigar_buf = Calloc(bd->cigar_buf_sz, char);
    return bd;
}

BAM_DATA
_init_BAM_DATA(SEXP ext, SEXP regions, SEXP flag, SEXP isSimpleCigar,
               SEXP tagFilter, SEXP mapqFilter, int reverseComplement,
               int yieldSize, int obeyQname, int asMates,
               char qnamePrefixEnd, char qnameSuffixStart, void *extra)
{
    int nrange = R_NilValue == regions ? 1 : LENGTH(VECTOR_ELT(regions, 0));
    BAM_DATA bd =
        _Calloc_BAM_DATA(1 == nrange ?
                         5 * BAM_INIT_SIZE : BAM_INIT_SIZE, 32768);
    bd->parse_status = BAM_PARSE_STATUS_OK;
    bd->bfile = BAMFILE(ext);
    bd->irange = bd->bfile->irange0;
    bd->nrange = nrange;
    bd->irec = bd->iparsed = 0;
    bd->keep_flag[0] = INTEGER(flag)[0];
    bd->keep_flag[1] = INTEGER(flag)[1];
    bd->cigar_flag = LOGICAL(isSimpleCigar)[0];
    bd->tagfilter = _tagFilter_as_C_types(tagFilter);
    int mapqfilter = INTEGER(mapqFilter)[0];
    bd->mapqfilter = mapqfilter == NA_INTEGER ? 0 : mapqfilter; /* uint32_t */
    bd->reverseComplement = reverseComplement;
    bd->yieldSize = yieldSize;
    bd->obeyQname = obeyQname;
    bd->asMates = asMates;
    bd->qnamePrefixEnd = qnamePrefixEnd;
    bd->qnameSuffixStart = qnameSuffixStart;
    bd->extra = extra;
    return bd;
}

void _Free_BAM_DATA(BAM_DATA bd)
{
    _Free_C_TAGFILTER(bd->tagfilter);
    Free(bd->cigar_buf);
    Free(bd);
}


static void _grow_BAM_DATA_cigar(BAM_DATA bd)
{
    bd->cigar_buf_sz = bd->cigar_buf_sz * 1.6;
    bd->cigar_buf = Realloc(bd->cigar_buf, bd->cigar_buf_sz, char);
}

BAM_FILE _bam_file_BAM_DATA(BAM_DATA bd)
{
    return bd->bfile;
}

/* count */

int _count1_BAM_DATA(const bam1_t * bam, BAM_DATA bd)
{
    bd->irec += 1;
    if (!_filter1_BAM_DATA(bam, bd))
        return 0;
    SEXP cnt = (SEXP) (bd->extra);
    INTEGER(VECTOR_ELT(cnt, 0))[bd->irange] += 1;
    REAL(VECTOR_ELT(cnt, 1))[bd->irange] += bam->core.l_qseq;
    bd->iparsed += 1;
    return 1;
}

/* parse helpers */

static char *_bamseq(const bam1_t * bam, BAM_DATA bd)
{
    static const char key[] = {
        '-', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
        'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'
    };

    const uint32_t len = bam->core.l_qseq;
    const unsigned char *seq = bam1_seq(bam);
    char *s = Calloc(len + 1, char);
    for (uint32_t i = 0; i < len; ++i)
        s[i] = key[bam1_seqi(seq, i)];
    if (bd->reverseComplement && (bam1_strand(bam) == 1))
        _reverseComplement(s, len);
    s[len] = '\0';
    return s;
}

static char *_bamqual(const bam1_t * bam, BAM_DATA bd)
{
    const uint32_t len = bam->core.l_qseq;
    const unsigned char *bamq = bam1_qual(bam);
    char *s = Calloc(len + 1, char);
    for (uint32_t i = 0; i < len; ++i)
        s[i] = bamq[i] + 33;
    if (bd->reverseComplement && (bam1_strand(bam) == 1))
        _reverse(s, len);
    s[len] = '\0';
    return s;
}

static const char *_map(khash_t(str) * h, const char *s)
{
    khiter_t k = kh_get(str, h, s);
    if (kh_end(h) == k) {
        int ret;
        char *buf = Calloc(strlen(s) + 1, char);
        if (!buf)
            Rf_error("_map: failed to allocate memory");
        strcpy(buf, s);
        k = kh_put(str, h, buf, &ret);
    }
    return kh_key(h, k);
}

static int _bamcigar(const uint32_t * cigar, const uint32_t n_cigar, char *buf,
                     int buf_sz)
{
    const char lookup[] = BAM_CIGAR_STR;
    buf[0] = '\0';
    for (uint32_t i = 0; i < n_cigar; ++i) {
        int n = snprintf(buf, buf_sz, "%u%c", cigar[i] >> 4,
                         lookup[cigar[i] & BAM_CIGAR_MASK]);
        if (n >= buf_sz)
            return -1;
        buf += n;
        buf_sz -= n;
    }
    return buf_sz;
}

static void _tag_type_check(const char *tagname, SEXP tag, SEXPTYPE is)
{
    SEXPTYPE was = TYPEOF(tag);
    if (was == is)
        return;
    error("tag '%s' type is inconsistent; was '%s', is '%s'",
          tagname, Rf_type2char(was), Rf_type2char(is));
}

static SEXP _bamtags_B(uint8_t *s) {
    uint8_t sub_type = *(s++);
    int32_t n;
    memcpy(&n, s, 4);

    SEXP tag_B = R_NilValue;
    switch(sub_type) {
    case 'c':
    case 'C':
    case 'i':
    case 'I':
    case 's':
    case 'S':
        tag_B = NEW_INTEGER(n); PROTECT(tag_B);
        for (int j = 0; j < n; ++j)
            INTEGER(tag_B)[j] = NA_INTEGER;
        break;
    case 'f':
        tag_B = NEW_NUMERIC(n); PROTECT(tag_B);
        for (int j = 0; j < n; ++j)
            REAL(tag_B)[j] = NA_REAL;
        break;
    }

    s += 4; // now point to the start of the array
    for (int i = 0; i < n; ++i) { // FIXME: better performance: loop after "if"
        if ('c' == sub_type || 'C' == sub_type) {
            INTEGER(tag_B)[i] = (int) (*(uint8_t *) s);
            ++s;
        } else if ('s' == sub_type || 'S' == sub_type) {
            INTEGER(tag_B)[i] = (int) (*(int16_t *) s);
            s += 2;
        } else if ('i' == sub_type || 'I' == sub_type) {
            INTEGER(tag_B)[i] = (int) (*(int32_t *) s);
            s += 4;
        } else if ('f' == sub_type) {
            REAL(tag_B)[i] = *(float *) s;
            s += 4;
        }
    }

    UNPROTECT(1);
    return tag_B;
}

static void _bamtags(const bam1_t * bam, BAM_DATA bd, SEXP tags)
{
    SCAN_BAM_DATA sbd = (SCAN_BAM_DATA) bd->extra;
    static char buf_A[2];
    int idx = sbd->icnt;
    SEXP nms = GET_ATTR(tags, R_NamesSymbol);
    buf_A[1] = '\0';            /* strictly necessary only once */

    for (int i = 0; i < LENGTH(nms); ++i) {
        const char *tagname = CHAR(STRING_ELT(nms, i));
        uint8_t *aux = bam_aux_get(bam, tagname);
        if (0 == aux)
            continue;           /* no matching tag found */
        SEXP tag = VECTOR_ELT(tags, i);
        if (R_NilValue == tag) {	/* allocate */
            int n = sbd->ncnt;
            switch (aux[0]) {
            case 'c':
            case 'C':
            case 'i':
            case 'I':
            case 's':
            case 'S':
                tag = NEW_INTEGER(n);
                for (int j = 0; j < n; ++j)
                    INTEGER(tag)[j] = NA_INTEGER;
                break;
            case 'd':
            case 'f':
                tag = NEW_NUMERIC(n);
                for (int j = 0; j < n; ++j)
                    REAL(tag)[j] = NA_REAL;
                break;
            case 'A':
            case 'Z':
                tag = NEW_CHARACTER(n);
                for (int j = 0; j < n; ++j)
                    SET_STRING_ELT(tag, j, NA_STRING);
                break;
            case 'H':
                tag = NEW_RAW(n);
                break;
            case 'B':
                tag = NEW_LIST(n);
                break;
            default:
                error("unknown tag type '%c'", aux[0]);
                break;
            }
            SET_VECTOR_ELT(tags, i, tag);
        }
        switch (aux[0]) {
        case 'c':
        case 'C':
        case 's':
        case 'S':
        case 'i':
        case 'I':
            _tag_type_check(tagname, tag, INTSXP);
            INTEGER(tag)[idx] = bam_aux2i(aux);
            break;
        case 'f':
        case 'd':
            _tag_type_check(tagname, tag, REALSXP);
            REAL(tag)[idx] = bam_aux2f(aux);
            break;
        case 'A':
            _tag_type_check(tagname, tag, STRSXP);
            sprintf(buf_A, "%c", bam_aux2A(aux));
            SET_STRING_ELT(tag, idx, mkChar(buf_A));
            break;
        case 'Z':
            _tag_type_check(tagname, tag, STRSXP);
            SET_STRING_ELT(tag, idx, mkChar(bam_aux2Z(aux)));
            break;
        case 'H':              /* FIXME: one byte or many? */
            _tag_type_check(tagname, tag, RAWSXP);
            RAW(tag)[idx] = aux[1];
            break;
        case 'B':
            _tag_type_check(tagname, tag, VECSXP);
            SET_VECTOR_ELT(tag, idx, _bamtags_B(&aux[1]));
            break;
        default:
            error("unknown tag type '%c'", aux[0]);
            break;
        }
    }
}

/* parse */

int _filter1_BAM_DATA(const bam1_t * bam, BAM_DATA bd)
{
    
    /* tagfilter */
    if (bd->tagfilter != NULL && !_tagfilter(bam, bd->tagfilter, bd->irec))
        return 0;
    if (bam->core.qual < bd->mapqfilter)
        return 0;

    /*
       flag : 1101
       keep0: 1111
       keep1: 1111
       test = (keep0 & ~flag) | (keep1 & flag) = 0010 | 1101 = 1111
       ~test = 0000 = FALSE

       flag:  1101
       keep0: 1101
       keep1: 1111
       test = (keep0 & ~flag) | (keep1 & flag) = 0000 | 1101 = 1101
       ~test = 0010 = TRUE

       flag:  1101
       keep0: 1111
       keep1: 1101
       test = (keep0 & ~flag) | (keep1 & flag) = 0010 | 1101 = 1111
       ~test = 0000 = FALSE
     */

    uint32_t test = (bd->keep_flag[0] & ~bam->core.flag) |
        (bd->keep_flag[1] & bam->core.flag);
    if (~test & 2047u)
        return 0;

    uint32_t *cigar = bam1_cigar(bam);
    uint32_t n_cigar = bam->core.n_cigar;
    if (bd->cigar_flag == CIGAR_SIMPLE) {
        if (!(n_cigar == 0 ||
              (n_cigar == 1 && ((cigar[0] & BAM_CIGAR_MASK) == 0))))
            return 0;
    }
    return 1;
}

int _filter_and_parse1_BAM_DATA(const bam1_t *bam, BAM_DATA bd)
{
    bd->irec += 1;
    if (!_filter1_BAM_DATA(bam, bd))
        return 0;
    return _parse1_BAM_DATA(bam, bd);
}

int _parse1_BAM_DATA(const bam1_t *bam, BAM_DATA bd)
{
    SCAN_BAM_DATA sbd = (SCAN_BAM_DATA) bd->extra;
    SEXP r = _get_or_grow_SCAN_BAM_DATA(bd, -1), s;
    int idx = sbd->icnt;
    char *buf;

    for (int i = 0; i < LENGTH(r); ++i) {
        if (R_NilValue == (s = VECTOR_ELT(r, i)))
            continue;
        switch (i) {
        case QNAME_IDX:
            buf = Calloc(strlen(bam1_qname(bam)) + 1, char);
            if (!buf)
                Rf_error("_parse1: failed to allocate memory");
            strcpy(buf, bam1_qname(bam));
            sbd->qname[idx] = buf;
            break;
        case FLAG_IDX:
            sbd->flag[idx] = bam->core.flag;
            break;
        case RNAME_IDX:
            sbd->rname[idx] =
                bam->core.tid < 0 ? NA_INTEGER : bam->core.tid + 1;
            break;
        case STRAND_IDX:
            sbd->strand[idx] = bam->core.flag & BAM_FUNMAP ?
                NA_INTEGER : (bam1_strand(bam) + 1);
            break;
        case POS_IDX:
            sbd->pos[idx] = bam->core.flag & BAM_FUNMAP ?
                NA_INTEGER : bam->core.pos + 1;
            break;
        case QWIDTH_IDX:
            sbd->qwidth[idx] = bam->core.flag & BAM_FUNMAP ?
                NA_INTEGER : bam_cigar2qlen(&bam->core, bam1_cigar(bam));
            break;
        case MAPQ_IDX:
            if ((bam->core.flag & BAM_FUNMAP))
                sbd->mapq[idx] = NA_INTEGER;
            else sbd->mapq[idx] = bam->core.qual;
            break;
        case CIGAR_IDX:
            if (bam->core.flag & BAM_FUNMAP)
                sbd->cigar[idx] = NULL;
            else {
                while (_bamcigar(bam1_cigar(bam), bam->core.n_cigar,
                              bd->cigar_buf, bd->cigar_buf_sz) < 0)
                    _grow_BAM_DATA_cigar(bd);
                sbd->cigar[idx] = _map(sbd->cigarhash, bd->cigar_buf);
            }
            break;
        case MRNM_IDX:
            sbd->mrnm[idx] = bam->core.mtid < 0 ?
                NA_INTEGER : bam->core.mtid + 1;
            break;
        case MPOS_IDX:
            sbd->mpos[idx] = bam->core.flag & BAM_FMUNMAP ?
                NA_INTEGER : bam->core.mpos + 1;
            break;
        case ISIZE_IDX:
            sbd->isize[idx] =
                bam->core.flag & (BAM_FUNMAP | BAM_FMUNMAP) ?
                NA_INTEGER : bam->core.isize;
            break;
        case SEQ_IDX:
            sbd->seq[idx] = _bamseq(bam, bd);
            break;
        case QUAL_IDX:
            sbd->qual[idx] = _bamqual(bam, bd);
            break;
        case TAG_IDX:
            _bamtags(bam, bd, s);
            break;
        case PARTITION_IDX:
            sbd->partition[idx] = sbd->partition_id;
            break;
        case MATES_IDX:
            sbd->mates[idx] = sbd->mates_flag;
            break;
        default:
            Rf_error("[Rsamtools internal]: unhandled _parse1");
            break;
        }
    }
    sbd->icnt += 1;
    bd->iparsed += 1;
    return 1;
}

void _finish1range_BAM_DATA(BAM_DATA  bd)
{
    bam_header_t *header = _bam_file_BAM_DATA(bd)->file->header;
    SCAN_BAM_DATA sbd = (SCAN_BAM_DATA) bd->extra;
    _finish1range_SCAN_BAM_DATA(sbd, header, bd->irange);
}

