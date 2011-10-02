#include "samtools/sam.h"
#include "io_sam.h"
#include "bamfile.h"
#include "encode.h"
#include "utilities.h"
#include "IRanges_interface.h"
#include "Biostrings_interface.h"

/* from samtoools/bam_sort.c */
void
bam_sort_core(int is_by_qname, const char *fn, const char *prefix,
              size_t max_mem);

typedef enum {
    OK = 0, SEQUENCE_BUFFER_ALLOCATION_ERROR=1,
    CIGAR_BUFFER_OVERFLOW_ERROR = 2
} _BAM_PARSE_STATUS;

typedef struct {
    int BLOCKSIZE;            /* size to grow vectors */
    char *CIGAR_BUF;          /* string representation of CIGAR */
    uint32_t CIGAR_BUF_SZ;

    _BAM_PARSE_STATUS parse_status;
    _BAM_FILE *bfile;
    int idx, icnt, ncnt, irange, nrange;
    uint32_t keep_flag[2], cigar_flag;
    Rboolean reverseComplement;
    const char **seq, **qual;

    void *extra;
} _BAM_DATA;

typedef int (*_PARSE1_FUNC)(const bam1_t *, void *);
typedef void (_FINISH1_FUNC)(_BAM_DATA *);

static const char *TMPL_ELT_NMS[] = {
    "qname", "flag", "rname", "strand", "pos", "qwidth", "mapq", "cigar",
    "mrnm", "mpos", "isize", "seq", "qual", "tag"
    /* "vtype", "value" */
};

static const int N_TMPL_ELTS = sizeof(TMPL_ELT_NMS) / sizeof(const char *);

enum {
    QNAME_IDX = 0, FLAG_IDX, RNAME_IDX, STRAND_IDX, POS_IDX, QWIDTH_IDX,
    MAPQ_IDX, CIGAR_IDX, MRNM_IDX, MPOS_IDX, ISIZE_IDX, SEQ_IDX,
    QUAL_IDX, TAG_IDX
};

enum { CIGAR_SIMPLE = 1 };

SEXP _count_bam(SEXP bfile, SEXP space, SEXP keepFlags,
		SEXP isSimpleCigar);

void _scan_bam_get1range_tag(SEXP tags, int len);

SEXP _scan_bam_get1range(_BAM_DATA *bd, int len);

void
_bam_check_template_list(SEXP template_list)
{
    if (!IS_LIST(template_list) ||
        LENGTH(template_list) != N_TMPL_ELTS)
        Rf_error("'template' must be list(%d)", N_TMPL_ELTS);
    SEXP names = GET_ATTR(template_list, R_NamesSymbol);
    if (!IS_CHARACTER(names) || LENGTH(names) != N_TMPL_ELTS)
        Rf_error("'names(template)' must be character(%d)",
                 N_TMPL_ELTS);
    for (int i = 0; i < LENGTH(names); ++i)
        if (strcmp(TMPL_ELT_NMS[i], CHAR(STRING_ELT(names, i))) != 0)
            Rf_error("'template' names do not match scan_bam_template\n'");
}

_BAM_DATA *
_Calloc_BAM_DATA(int blocksize, int cigar_buf_sz)
{
    _BAM_DATA *bd = (_BAM_DATA *) Calloc(1, _BAM_DATA);
    bd->BLOCKSIZE = blocksize;
    bd->CIGAR_BUF_SZ = cigar_buf_sz;
    bd->CIGAR_BUF = Calloc(bd->CIGAR_BUF_SZ, char);
    bd->seq = bd->qual = NULL;

    return bd;
}

void
_Free_BAM_DATA(_BAM_DATA *bd)
{
    Free(bd->CIGAR_BUF);
    Free(bd->seq);
    Free(bd->qual);
    Free(bd);
}

static char *
_bamseq(const bam1_t *bam, _BAM_DATA *bd)
{
    static const char key[] = {
        '\0', 'A', 'C',  '\0',  'G', '\0', '\0', '\0',
        'T', '\0', '\0', '\0', '\0', '\0', '\0', 'N'
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

static char *
_bamqual(const bam1_t *bam, _BAM_DATA *bd)
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

int
_bamcigar(const uint32_t *cigar, const uint32_t n_cigar,
          char *buf, int buf_sz)
{
    const char lookup[] = { 'M', 'I', 'D', 'N', 'S', 'H', 'P' };
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

void
_tag_type_check(const char *tagname, SEXP tag, SEXPTYPE is)
{
    SEXPTYPE was = TYPEOF(tag);
    if (was == is)
        return;
    error("tag '%s' type is inconsistent; was '%s', is '%s'",
          tagname, Rf_type2char(was), Rf_type2char(is));
}

void
_bamtags(const bam1_t *bam, _BAM_DATA *bd, SEXP tags)
{
    static char *buf_A;
    int idx = bd->icnt;
    SEXP nms = GET_ATTR(tags, R_NamesSymbol);
    for (int i = 0; i < LENGTH(nms); ++i) {
        const char *tagname = CHAR(STRING_ELT(nms, i));
        uint8_t *aux = bam_aux_get(bam, tagname);
        if (0 == aux)
            continue;           /* no matching tag found */
        SEXP tag = VECTOR_ELT(tags, i);
        if (R_NilValue == tag) { /* allocate */
            int n = bd->ncnt;
            switch (aux[0]) {
            case 'c': case 'C': case 'i': case 'I': case 's': case 'S':
                tag = NEW_INTEGER(n);
                for (int j = 0; j < n; ++j)
                    INTEGER(tag)[j] = NA_INTEGER;
                break;
            case 'd': case 'f':
                tag = NEW_NUMERIC(n);
                for (int j = 0; j < n; ++j)
                    REAL(tag)[j] = NA_REAL;
                break;
            case 'A': case 'Z':
                tag = NEW_CHARACTER(n);
                for (int j = 0; j < n; ++j)
                    SET_STRING_ELT(tag, j, NA_STRING);
                if ('A' == aux[0]) {
                    buf_A = R_alloc(2, sizeof(char));
                    buf_A[1] = '\0';
                }
                break;
            case 'H':
                tag = NEW_RAW(n);
                break;
            default:
                error("unknown tag type '%c'", aux[0]);
                break;
            }
            SET_VECTOR_ELT(tags, i, tag);
        }
        switch (aux[0]) {
        case 'c': case 'C': case 's': case 'S': case 'i': case 'I':
            _tag_type_check(tagname, tag, INTSXP);
            INTEGER(tag)[idx] = bam_aux2i(aux);
            break;
        case 'f':
            _tag_type_check(tagname, tag, REALSXP);
            REAL(tag)[idx] = (double) bam_aux2f(aux);
            break;
        case 'd':
            _tag_type_check(tagname, tag, REALSXP);
            REAL(tag)[idx] = bam_aux2d(aux);
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
        case 'H':               /* FIXME: one byte or many? */
            _tag_type_check(tagname, tag, RAWSXP);
            RAW(tag)[idx] = aux[1];
            break;
        default:
            error("unknown tag type '%c'", aux[0]);
            break;
        }
    }
}

/* parse header */

SEXP
_read_bam_header(SEXP ext)
{
    samfile_t *sfile = BAMFILE(ext)->file;
    bam_header_t *header = sfile->header;
    int n_elts = header->n_targets;

    SEXP ans = PROTECT(NEW_LIST(2));

    /* target length / name */
    SET_VECTOR_ELT(ans, 0, NEW_INTEGER(n_elts));
    SEXP tlen = VECTOR_ELT(ans, 0); /* target length */
    SEXP tnm = PROTECT(NEW_CHARACTER(n_elts)); /* target name */
    setAttrib(tlen, R_NamesSymbol, tnm);
    UNPROTECT(1);
    for (int j = 0; j < n_elts; ++j) {
        INTEGER(tlen)[j] = header->target_len[j];
        SET_STRING_ELT(tnm, j, mkChar(header->target_name[j]));
    }

    /* 'aux' character string */
    char *txt = (char *) R_alloc(header->l_text + 1, sizeof(char));
    strncpy(txt, header->text, header->l_text);
    txt[header->l_text] = '\0';
    SET_VECTOR_ELT(ans, 1, mkString(txt));

    SEXP nms = PROTECT(NEW_CHARACTER(2));
    SET_STRING_ELT(nms, 0, mkChar("targets"));
    SET_STRING_ELT(nms, 1, mkChar("text"));
    setAttrib(ans, R_NamesSymbol, nms);
    UNPROTECT(2);
    return ans;
}

_BAM_DATA *
_init_BAM_DATA(SEXP ext, SEXP space, SEXP flag, SEXP isSimpleCigar,
               Rboolean reverseComplement)
{
    int nrange =
        R_NilValue == space ? 1 : LENGTH(VECTOR_ELT(space, 0));
    _BAM_DATA *bd;
    bd = _Calloc_BAM_DATA(1 == nrange ? 5 * 1048576 : 1048576, 32768);
    bd->parse_status = 0;
    bd->bfile = BAMFILE(ext);
    bd->irange = 0;
    bd->nrange = nrange;
    bd->idx = 0;
    bd->keep_flag[0] = INTEGER(flag)[0];
    bd->keep_flag[1] = INTEGER(flag)[1];
    bd->cigar_flag = LOGICAL(isSimpleCigar)[0];
    bd->reverseComplement = reverseComplement;

    bd->icnt = bd->ncnt = 0;
    bd->seq = bd->qual = NULL;
    return bd;
}

Rboolean
_bam_filter(const bam1_t *bam, _BAM_DATA *bd)
{
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
        return FALSE;

    uint32_t *cigar = bam1_cigar(bam);
    uint32_t n_cigar = bam->core.n_cigar;
    if (bd->cigar_flag == CIGAR_SIMPLE)
    {
        if (!(n_cigar == 0 ||
              (n_cigar == 1  && ((cigar[0] & BAM_CIGAR_MASK) == 0))))
            return FALSE;
    }
    return TRUE;
}

int
_scan_bam_all(_BAM_DATA *bd, _PARSE1_FUNC parse1, _FINISH1_FUNC finish1)
{
    bam1_t *bam = bam_init1();
    int r = 0;

    bam_seek(bd->bfile->file->x.bam, bd->bfile->pos0, SEEK_SET);
    while ((r = samread(bd->bfile->file, bam)) >= 0) {
        int result = (*parse1)(bam, bd);
        if (result < 0) return result;
    }
    if (NULL != finish1)
        (*finish1)(bd);
    return bd->icnt;
}

int
_scan_bam_fetch(_BAM_DATA *bd, SEXP space, int* start, int* end,
                _PARSE1_FUNC parse1, _FINISH1_FUNC finish1)
{
    int tid;
    samfile_t *sfile = bd->bfile->file;
    bam_index_t *bindex = bd->bfile->index;
    int n_tot = 0;

    for (int irange = 0; irange < LENGTH(space); ++irange) {
        const char* spc = translateChar(STRING_ELT(space, irange));
        const int starti =
            start[irange] > 0 ? start[irange] - 1 : start[irange];
        for (tid = 0; tid < sfile->header->n_targets; ++tid) {
            if (strcmp(spc, sfile->header->target_name[tid]) == 0)
                break;
        }
        if (tid == sfile->header->n_targets) {
            Rf_warning("space '%s' not in BAM header", spc);
            return -1;
        }

        bam_fetch(sfile->x.bam, bindex, tid,
                  starti, end[irange], bd, parse1);
        n_tot += bd->icnt;
        if (NULL != finish1)
            (*finish1)(bd);
        bd->irange += 1;
    }

    return n_tot;
}

int
_do_scan_bam(_BAM_DATA *bd, SEXP space, _PARSE1_FUNC parse1,
             _FINISH1_FUNC finish1)
{
    int status;

    if (R_NilValue == space)
	/* everything */
        status = _scan_bam_all(bd, parse1, finish1);
    else {                   /* fetch */
	if (NULL == bd->bfile->index)
	    Rf_error("valid 'index' file required");
        status = _scan_bam_fetch(bd, VECTOR_ELT(space, 0),
                                 INTEGER(VECTOR_ELT(space, 1)),
                                 INTEGER(VECTOR_ELT(space, 2)),
                                 parse1, finish1);
    }

    return status;
}

/* parse */

static SEXP
_get_lkup(const char *baseclass)
{
    SEXP lkup = R_NilValue;
    switch(*baseclass) {
    case 'D':
        lkup = _get_encoding_lookup("B", "DNA");
        break;
    case 'B':
        break;
    default:
        Rf_error("Rsamtools internal: '%s' unhandled in _get_lkup",
                 baseclass);
        break;
    }
    return lkup;
}

static SEXP
_tmpl_DNAStringSet()
{
    CharAEAE aeae = new_CharAEAE(0, 0);
    SEXP lkup = PROTECT(_get_lkup("DNAString"));
    SEXP ans = new_XRawList_from_CharAEAE("DNAStringSet", "DNAString",
                                          &aeae, lkup);
    UNPROTECT(1);
    return ans;
}

static SEXP
_tmpl_BStringSet()
{
    CharAEAE aeae = new_CharAEAE(0, 0);
    return new_XRawList_from_CharAEAE("BStringSet", "BString", &aeae,
                                      R_NilValue);
}

static SEXP
_tmpl_PhredQuality()
{
    SEXP xstringset, s, t, nmspc, result;
    PROTECT(xstringset = _tmpl_BStringSet());
    PROTECT(nmspc = _get_namespace("Rsamtools"));
    NEW_CALL(s, t, "PhredQuality", nmspc, 2);
    CSET_CDR(t, "x", xstringset);
    CEVAL_TO(s, nmspc, result);
    UNPROTECT(2);
    return result;
}

SEXP
scan_bam_template(SEXP tag)
{
    if (R_NilValue != tag)
        if (!IS_CHARACTER(tag))
            Rf_error("'tag' must be NULL or 'character()'");
    SEXP tmpl = PROTECT(NEW_LIST(N_TMPL_ELTS));
    SET_VECTOR_ELT(tmpl, QNAME_IDX, NEW_CHARACTER(0));
    SET_VECTOR_ELT(tmpl, FLAG_IDX, NEW_INTEGER(0));
    SET_VECTOR_ELT(tmpl, RNAME_IDX, NEW_INTEGER(0));
    SET_VECTOR_ELT(tmpl, STRAND_IDX, NEW_INTEGER(0));
    SET_VECTOR_ELT(tmpl, POS_IDX, NEW_INTEGER(0));
    SET_VECTOR_ELT(tmpl, QWIDTH_IDX, NEW_INTEGER(0));
    SET_VECTOR_ELT(tmpl, MAPQ_IDX, NEW_INTEGER(0));
    SET_VECTOR_ELT(tmpl, CIGAR_IDX, NEW_CHARACTER(0));
    SET_VECTOR_ELT(tmpl, MRNM_IDX, NEW_INTEGER(0));
    SET_VECTOR_ELT(tmpl, MPOS_IDX, NEW_INTEGER(0));
    SET_VECTOR_ELT(tmpl, ISIZE_IDX, NEW_INTEGER(0));
    SET_VECTOR_ELT(tmpl, SEQ_IDX, _tmpl_DNAStringSet());
    SET_VECTOR_ELT(tmpl, QUAL_IDX, _tmpl_PhredQuality());
    if (R_NilValue == tag) {
        SET_VECTOR_ELT(tmpl, TAG_IDX, R_NilValue);
    } else {
        SET_VECTOR_ELT(tmpl, TAG_IDX, NEW_LIST(LENGTH(tag)));
        SET_ATTR(VECTOR_ELT(tmpl, TAG_IDX), R_NamesSymbol, tag);
    }

    SEXP names = PROTECT(NEW_CHARACTER(N_TMPL_ELTS));
    for (int i = 0; i < N_TMPL_ELTS; ++i)
        SET_STRING_ELT(names, i, mkChar(TMPL_ELT_NMS[i]));
    SET_ATTR(tmpl, R_NamesSymbol, names);
    UNPROTECT(2);
    return tmpl;
}

int
_scan_bam_parse1(const bam1_t *bam, void *data)
{
    _BAM_DATA *bd = (_BAM_DATA *) data;
    bd->idx += 1;
    if (FALSE == _bam_filter(bam, bd))
        return 0;

    SEXP r = _scan_bam_get1range(bd, -1), s;
    int idx = bd->icnt;

    for (int i = 0; i < LENGTH(r); ++i) {
        if ((s = VECTOR_ELT(r, i)) == R_NilValue)
            continue;
        switch(i) {
        case QNAME_IDX:
            SET_STRING_ELT(s, idx, mkChar(bam1_qname(bam)));
            break;
        case FLAG_IDX:
            INTEGER(s)[idx] = bam->core.flag;
            break;
        case RNAME_IDX:
            INTEGER(s)[idx] =
                bam->core.tid < 0 ? NA_INTEGER : bam->core.tid + 1;
            break;
        case STRAND_IDX:
            INTEGER(s)[idx] = bam->core.flag & BAM_FUNMAP ?
                NA_INTEGER : (bam1_strand(bam) + 1);
            break;
        case POS_IDX:
            INTEGER(s)[idx] = bam->core.flag & BAM_FUNMAP ?
                NA_INTEGER : bam->core.pos + 1;
            break;
        case QWIDTH_IDX:
            INTEGER(s)[idx] = bam->core.flag & BAM_FUNMAP ?
                NA_INTEGER : bam_cigar2qlen(&bam->core, bam1_cigar(bam));
            break;
        case MAPQ_IDX:
            INTEGER(s)[idx] = bam->core.flag & BAM_FUNMAP ?
                NA_INTEGER : bam->core.qual;
            break;
        case CIGAR_IDX:
            if (bam->core.flag & BAM_FUNMAP)
                SET_STRING_ELT(s, idx, NA_STRING);
            else {
                if (_bamcigar(bam1_cigar(bam), bam->core.n_cigar,
                              bd->CIGAR_BUF, bd->CIGAR_BUF_SZ) < 0)
                {
                    bd->parse_status |= CIGAR_BUFFER_OVERFLOW_ERROR;
                    return -bd->icnt;
                }
                SET_STRING_ELT(s, idx, mkChar(bd->CIGAR_BUF));
            }
            break;
        case MRNM_IDX:
            INTEGER(s)[idx] = bam->core.mtid < 0 ?
                NA_INTEGER : bam->core.mtid + 1;
            break;
        case MPOS_IDX:
            INTEGER(s)[idx] = bam->core.flag & BAM_FMUNMAP ?
                NA_INTEGER : bam->core.mpos + 1;
            break;
        case ISIZE_IDX:
            INTEGER(s)[idx] =
                bam->core.flag & (BAM_FUNMAP | BAM_FMUNMAP) ?
                NA_INTEGER : bam->core.isize;
            break;
        case SEQ_IDX:
            bd->seq[idx] = _bamseq(bam, bd);
            break;
        case QUAL_IDX:
            bd->qual[idx] = _bamqual(bam, bd);
            break;
        case TAG_IDX:
            _bamtags(bam, bd, s);
            break;
        default:
            break;
        }
    }
    bd->icnt += 1;
    return 1;
}

SEXP
_as_XStringSet(const char **key, int len, const char *baseclass)
{
    char classname[40];  /* longest string should be "DNAStringSet" */

    if (snprintf(classname, sizeof(classname), "%sSet", baseclass)
        >= (int) sizeof(classname))
        error("Rsamtools internal error in _as_XStringSet(): "
              "'classname' buffer too small");

    SEXP lkup = _get_lkup(baseclass);
    int *lkup0, lkup_length = 0;
    if (R_NilValue == lkup) {
        lkup0 = NULL;
    } else {
        lkup0 = INTEGER(lkup);
        lkup_length = LENGTH(lkup);
    }

    SEXP width, ans;
    int i;
    PROTECT(width = NEW_INTEGER(len));
    for (i = 0; i < len; ++i)
        INTEGER(width)[i] = strlen(key[i]);
    PROTECT(ans = alloc_XRawList(classname, baseclass, width));

    cachedXVectorList cache;
    cachedCharSeq dest;
    cache = cache_XVectorList(ans);
    for (i = 0; i < len; ++i) {
        const char *seq = key[i];
        dest = get_cachedXRawList_elt(&cache, i);
        Ocopy_bytes_to_i1i2_with_lkup(
            0, dest.length - 1, (char *) dest.seq, dest.length,
            seq, strlen(seq), lkup0, lkup_length);
    }

    UNPROTECT(2);
    return ans;
}

SEXP
_as_PhredQuality(const char **key, int len)
{
    SEXP xstringset = PROTECT(_as_XStringSet(key, len, "BString"));

    SEXP s, t, nmspc, result;
    nmspc = PROTECT(_get_namespace("Rsamtools"));
    NEW_CALL(s, t, "PhredQuality", nmspc, 2);
    CSET_CDR(t, "x", xstringset);
    CEVAL_TO(s, nmspc, result);
    UNPROTECT(2);
    return result;
}

SEXP
_scan_bam_result_init(SEXP template_list, SEXP names, SEXP space)
{
    const int nrange =
        R_NilValue == space ? 1 : Rf_length(VECTOR_ELT(space, 0));
    int i;

    SEXP result = PROTECT(NEW_LIST(nrange));
    /* result:
       range1: tmpl1, tmpl2...
       range2: tmpl1, tmpl2...
       ...
    */
    for (int irange = 0; irange < nrange; ++irange)
    {
        SEXP tag = VECTOR_ELT(template_list, TAG_IDX);
        SEXP tmpl;
        if (R_NilValue == tag)
            tmpl = PROTECT(scan_bam_template(R_NilValue));
        else {
            SEXP nms = getAttrib(tag, R_NamesSymbol);
            tmpl = PROTECT(scan_bam_template(nms));
        }
        for (i = 0; i < LENGTH(names); ++i) {
            if (TAG_IDX == i)
                continue;
            else if (R_NilValue == VECTOR_ELT(template_list, i))
                SET_VECTOR_ELT(tmpl, i, R_NilValue);
        }
        SET_VECTOR_ELT(result, irange, tmpl);
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return result;
}

void
_scan_bam_get1range_tag(SEXP tags, int len)
{
    int i, j;
    for (i = 0; i < LENGTH(tags); ++i) {
        SEXP elt = VECTOR_ELT(tags, i);
        int len0 = LENGTH(elt);
        SET_VECTOR_ELT(tags, i, Rf_lengthgets(elt, len));
        switch(TYPEOF(elt)) {
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

SEXP
_scan_bam_get1range(_BAM_DATA *bd, int len)
{
    SEXP r = VECTOR_ELT((SEXP) bd->extra, bd->irange);
    if (len < 0) {
        if (bd->icnt < bd->ncnt)
            return r;
        len = bd->ncnt + bd->BLOCKSIZE;
    }

    SEXP s;
    for (int i = 0; i < LENGTH(r); ++i) {
        if (R_NilValue == (s = VECTOR_ELT(r, i)))
            continue;
        SEXP elt;
        switch(i) {
        case TAG_IDX:
            _scan_bam_get1range_tag(s, len);
            break;
        case SEQ_IDX:
            bd->seq = Realloc(bd->seq, len, const char *);
            break;
        case QUAL_IDX:
            bd->qual = Realloc(bd->qual, len, const char *);
            break;
        default:
            elt = Rf_lengthgets(s, len);
            SET_VECTOR_ELT(r, i, elt);
            break;
        }
    }

    bd->ncnt = len;
    return r;
}

void
_scan_bam_finish1range(_BAM_DATA *bd)
{
    int i;
    bam_header_t *header = bd->bfile->file->header;
    SEXP r, s;

    r = _scan_bam_get1range(bd, bd->icnt);

    if (R_NilValue != (s = VECTOR_ELT(r, STRAND_IDX))) {
        SEXP strand_lvls = PROTECT(_get_strand_levels());
        _as_factor_SEXP(s, strand_lvls);
        UNPROTECT(1);
    }

    if (R_NilValue != (s = VECTOR_ELT(r, RNAME_IDX))) {
        _as_factor(s, (const char **) header->target_name,
                   header->n_targets);
    }

    if (R_NilValue != (s = VECTOR_ELT(r, MRNM_IDX))) {
        _as_factor(s, (const char **) header->target_name,
                   header->n_targets);
    }

    if (R_NilValue != (s = VECTOR_ELT(r, SEQ_IDX)))
    {
        PROTECT(s = _as_XStringSet(bd->seq, bd->icnt, "DNAString"));
        SET_VECTOR_ELT(r, SEQ_IDX, s);
        UNPROTECT(1);

        for (i = 0; i < bd->icnt; ++i)
            Free(bd->seq[i]);
        Free(bd->seq);
        bd->seq = NULL;
    }

    if (R_NilValue != (s = VECTOR_ELT(r, QUAL_IDX)))
    {
        PROTECT(s = _as_PhredQuality(bd->qual, bd->icnt));
        SET_VECTOR_ELT(r, QUAL_IDX, s);
        UNPROTECT(1);

        for (i = 0; i < bd->icnt; ++i)
            Free(bd->qual[i]);
        Free(bd->qual);
        bd->qual = NULL;
    }

    bd->icnt = bd->ncnt = 0;
}

SEXP
_scan_bam(SEXP bfile, SEXP space, SEXP keepFlags, SEXP isSimpleCigar,
	  SEXP reverseComplement, SEXP template_list)
{
    SEXP names = PROTECT(GET_ATTR(template_list, R_NamesSymbol));
    _BAM_DATA *bd =
        _init_BAM_DATA(bfile, space, keepFlags, isSimpleCigar,
                       LOGICAL(reverseComplement)[0]);
    SEXP result =
        _scan_bam_result_init(template_list, names, space);
    PROTECT(result);
    bd->extra = (void *) result;

    int status = _do_scan_bam(bd, space, _scan_bam_parse1,
                              _scan_bam_finish1range);
    if (status < 0) {
        int idx = bd->idx;
        _BAM_PARSE_STATUS parse_status = bd->parse_status;
        _Free_BAM_DATA(bd);
        Rf_error("'scanBam' failed:\n  record: %d\n  error: %d",
                 idx, parse_status);
    }

    _Free_BAM_DATA(bd);
    UNPROTECT(2);
    return result;
}

/* count */

int
_count_bam1(const bam1_t *bam, void *data)
{
    _BAM_DATA *bd = (_BAM_DATA *) data;
    bd->idx += 1;
    if (FALSE == _bam_filter(bam, bd))
        return 0;
    SEXP cnt = (SEXP) (bd->extra);
    INTEGER(VECTOR_ELT(cnt, 0))[bd->irange] += 1;
    REAL(VECTOR_ELT(cnt, 1))[bd->irange] += bam->core.l_qseq;
    return 1;
}

SEXP
_count_bam(SEXP bfile, SEXP space, SEXP keepFlags, SEXP isSimpleCigar)
{
    _BAM_DATA *bd =
        _init_BAM_DATA(bfile, space, keepFlags, isSimpleCigar, FALSE);
    SEXP result = PROTECT(NEW_LIST(2));
    bd->extra = result;

    SET_VECTOR_ELT(result, 0, NEW_INTEGER(bd->nrange));
    SET_VECTOR_ELT(result, 1, NEW_NUMERIC(bd->nrange));
    for (int i = 0; i < bd->nrange; ++i) {
        INTEGER(VECTOR_ELT(result, 0))[i] =
            REAL(VECTOR_ELT(result, 1))[i] = 0;
    }

    SEXP nms = PROTECT(NEW_CHARACTER(2));
    SET_STRING_ELT(nms, 0, mkChar("records"));
    SET_STRING_ELT(nms, 1, mkChar("nucleotides"));
    setAttrib(result, R_NamesSymbol, nms);
    UNPROTECT(1);

    int status = _do_scan_bam(bd, space, _count_bam1, NULL);
    if (status < 0)
        result = R_NilValue;

    _Free_BAM_DATA(bd);
    UNPROTECT(1);
    return result;
}

void
scan_bam_cleanup()
{
    /* placeholder */
}

/* filterBam */

int
_filter_bam1(const bam1_t *bam, void *data)
{
    _BAM_DATA *bd = (_BAM_DATA *) data;
    bd->idx += 1;
    if (FALSE == _bam_filter(bam, bd))
        return 0;
    samwrite((samfile_t*) bd->extra, bam);
    return 1;
}

SEXP
_filter_bam(SEXP bfile, SEXP space, SEXP keepFlags,
	    SEXP isSimpleCigar, SEXP fout_name, SEXP fout_mode)
{
    /* open destination */
    _BAM_DATA *bd =
        _init_BAM_DATA(bfile, space, keepFlags, isSimpleCigar, FALSE);
    /* FIXME: this just copies the header... */
    bam_header_t *header = BAMFILE(bfile)->file->header;
    samfile_t *f_out =
        _bam_tryopen(translateChar(STRING_ELT(fout_name, 0)),
                     CHAR(STRING_ELT(fout_mode, 0)), header);
    bd->extra = f_out;

    int status = _do_scan_bam(bd, space, _filter_bam1, NULL);

    /* sort and index destintation ? */
    /* cleanup */
    samclose(f_out);
    _Free_BAM_DATA(bd);

    return status < 0 ? R_NilValue : fout_name;
}

/* sort_bam */

SEXP
sort_bam(SEXP filename, SEXP destination, SEXP isByQname,
	 SEXP maxMemory)
{
    if (!IS_CHARACTER(filename) || 1 != LENGTH(filename))
        Rf_error("'filename' must be character(1)");
    if (!IS_CHARACTER(destination) || 1 != LENGTH(destination))
        Rf_error("'destination' must be character(1)");
    if (!IS_LOGICAL(isByQname)  || LENGTH(isByQname) != 1)
        Rf_error("'isByQname' must be logical(1)");
    if (!IS_INTEGER(maxMemory) || LENGTH(maxMemory) != 1 ||
        INTEGER(maxMemory)[0] < 1)
        Rf_error("'maxMemory' must be a positive integer(1)");

    const char *fbam = translateChar(STRING_ELT(filename, 0));
    const char *fout = translateChar(STRING_ELT(destination, 0));
    int sortMode = asInteger(isByQname);

    size_t maxMem = (size_t) INTEGER(maxMemory)[0] * 1024 * 1024;
    bam_sort_core(sortMode, fbam, fout, maxMem);
    return destination;
}

/* index_bam */

SEXP
index_bam(SEXP indexname)
{
    if (!IS_CHARACTER(indexname) || 1 != LENGTH(indexname))
        Rf_error("'indexname' must be character(1)");
    const char *fbam = translateChar(STRING_ELT(indexname, 0));
    int status = bam_index_build(fbam);
    if (0 != status)
        Rf_error("failed to build index\n  file: %s", fbam);
    char *fidx = (char *) R_alloc(strlen(fbam) + 5, sizeof(char));
    sprintf(fidx, "%s.bai", fbam);
    return mkString(fidx);
}
