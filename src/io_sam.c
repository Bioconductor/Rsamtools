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
    int BUF_SZ, CIGAR_BUF_SZ; /* qual / seq and cigar scratch buffer */
    char *BUF, *CIGAR_BUF;    /* string representation of CIGAR */

    _BAM_PARSE_STATUS parse_status;
    _BAM_FILE *bfile;
    int nrec, idx, irange;
    uint32_t keep_flag[2], cigar_flag;
    Rboolean reverseComplement;

    CharAEAE *seq, *qual;

    void *extra;
} _BAM_DATA;

typedef int (*_PARSE1_FUNC)(const bam1_t *, void *);

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
_Calloc_BAM_DATA(int blocksize, int buf_sz, int cigar_buf_sz)
{
    _BAM_DATA *bd = (_BAM_DATA *) Calloc(1, _BAM_DATA);
    bd->BLOCKSIZE = blocksize;
    bd->BUF_SZ = buf_sz;
    bd->BUF = Calloc(bd->BUF_SZ, char);
    bd->CIGAR_BUF_SZ = cigar_buf_sz;
    bd->CIGAR_BUF = Calloc(bd->CIGAR_BUF_SZ, char);

    return bd;
}

void
_grow_BUF_BAM_DATA(_BAM_DATA *bd, int sz)
{
    Free(bd->BUF);
    bd->BUF_SZ = sz;
    bd->BUF = Calloc(bd->BUF_SZ, char);
}

void
_Free_BAM_DATA(_BAM_DATA *bd)
{
    Free(bd->BUF);
    Free(bd->CIGAR_BUF);
    Free(bd);
}

void
_bamseq(const bam1_t *bam, _BAM_DATA *bd, Rboolean reverseComplement)
{
    static const char key[] = {
        '\0', 'A', 'C',  '\0',  'G', '\0', '\0', '\0',
        'T', '\0', '\0', '\0', '\0', '\0', '\0', 'N'
    };

    uint32_t len = bam->core.l_qseq;
    unsigned char *seq = bam1_seq(bam);
    if (len + 1 > bd->BUF_SZ)
        _grow_BUF_BAM_DATA(bd, len + 1);
    for (int i = 0; i < len; ++i)
        bd->BUF[i] = key[bam1_seqi(seq, i)];
    if (reverseComplement && (bam1_strand(bam) == 1))
        _reverseComplement(bd->BUF, len);
    bd->BUF[len] = '\0';
    append_string_to_CharAEAE(bd->seq + bd->irange, bd->BUF);
}

void
_bamqual(const bam1_t *bam, _BAM_DATA *bd, Rboolean reverse)
{
    uint32_t len = bam->core.l_qseq;
    unsigned char *bamq = bam1_qual(bam);
    if (len + 1 > bd->BUF_SZ)
        _grow_BUF_BAM_DATA(bd, len + 1);
    for (int i = 0; i < len; ++i)
        bd->BUF[i] = bamq[i] + 33;
    if (reverse && (bam1_strand(bam) == 1))
        _reverse(bd->BUF, len);
    bd->BUF[len] = '\0';
    append_string_to_CharAEAE(bd->qual + bd->irange, bd->BUF);
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
    int idx = bd->idx;
    SEXP nms = GET_ATTR(tags, R_NamesSymbol);
    for (int i = 0; i < LENGTH(nms); ++i) {
        const char *tagname = CHAR(STRING_ELT(nms, i));
        uint8_t *aux = bam_aux_get(bam, tagname);
        if (0 == aux)
            continue;           /* no matching tag found */
        SEXP tag = VECTOR_ELT(tags, i);
        if (R_NilValue == tag) { /* allocate */
            int n = INTEGER(getAttrib(tags, install("count")))[0];
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
_init_BAM_DATA(SEXP ext, SEXP flag, SEXP isSimpleCigar,
               Rboolean reverseComplement)
{
    _BAM_DATA *bdata = _Calloc_BAM_DATA(1048576, 2048, 32768);
    bdata->parse_status = 0;
    bdata->bfile = BAMFILE(ext);
    bdata->nrec = bdata->idx = bdata->irange = 0;
    bdata->keep_flag[0] = INTEGER(flag)[0];
    bdata->keep_flag[1] = INTEGER(flag)[1];
    bdata->cigar_flag = LOGICAL(isSimpleCigar)[0];
    bdata->reverseComplement = reverseComplement;
    return bdata;
}

Rboolean
_bam_filter(const bam1_t *bam, _BAM_DATA *bdata)
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

    uint32_t test = (bdata->keep_flag[0] & ~bam->core.flag) |
        (bdata->keep_flag[1] & bam->core.flag);
    if (~test & 2047u)
        return FALSE;

    uint32_t *cigar = bam1_cigar(bam);
    uint32_t n_cigar = bam->core.n_cigar;
    if (bdata->cigar_flag == CIGAR_SIMPLE)
    {
        if (!(n_cigar == 0 ||
              (n_cigar == 1  && ((cigar[0] & BAM_CIGAR_MASK) == 0))))
            return FALSE;
    }
    return TRUE;
}

int
_scan_bam_all(_BAM_DATA *bd, _PARSE1_FUNC parse1)
{
    bam1_t *bam = bam_init1();
    int r = 0;
    while ((r = samread(bd->bfile->file, bam)) >= 0) {
        int result = (*parse1)(bam, bd);
        if (result < 0) return result;
    }
    return bd->idx;
}

int
_scan_bam_fetch(_BAM_DATA *bd, SEXP space, int* start, int* end,
                _PARSE1_FUNC parse1)
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
        n_tot += bd->idx;
        bd->irange += 1;
        bd->idx = 0;
    }

    return n_tot;
}

int
_do_scan_bam(_BAM_DATA *bdata, SEXP space, _PARSE1_FUNC parse1)
{
    int status;

    if (R_NilValue == space)	/* everything */
        status = _scan_bam_all(bdata, parse1);
    else {                   /* fetch */
	if (NULL == bdata->bfile->index)
	    Rf_error("valid 'index' file required");
        status = _scan_bam_fetch(bdata, VECTOR_ELT(space, 0),
				 INTEGER(VECTOR_ELT(space, 1)),
				 INTEGER(VECTOR_ELT(space, 2)), parse1);
    }

    return status;
}

/* parse */

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
    SET_VECTOR_ELT(tmpl, SEQ_IDX, NEW_CHARACTER(0));
    SET_VECTOR_ELT(tmpl, QUAL_IDX, NEW_CHARACTER(0));
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
    if (FALSE == _bam_filter(bam, bd))
        return 0;

    SEXP s, r = VECTOR_ELT((SEXP) bd->extra, bd->irange);
    int idx = bd->idx;
    Rboolean reverseComplement = bd->reverseComplement;

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
                    return -bd->idx;
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
            _bamseq(bam, bd, reverseComplement);
            break;
        case QUAL_IDX:
            _bamqual(bam, bd, reverseComplement);
            break;
        case TAG_IDX:
            _bamtags(bam, bd, s);
            break;
        default:
            break;
        }
    }
    bd->idx += 1;
    return 1;
}

SEXP
_scan_bam_result_init(SEXP count, SEXP template_list, SEXP names,
                      _BAM_DATA *bd)
{
    int nrange = LENGTH(VECTOR_ELT(count, 0)), i;
    SEXP result = PROTECT(NEW_LIST(nrange));
    /* result:
       range1: tmpl1, tmpl2...
       range2: tmpl1, tmpl2...
       ...
    */
    for (i = 0; i < LENGTH(names); ++i)
        if (SEQ_IDX == i)
            bd->seq = (CharAEAE *) R_alloc(nrange, sizeof(CharAEAE));
        else if (QUAL_IDX == i)
            bd->qual = (CharAEAE *) R_alloc(nrange, sizeof(CharAEAE));
    for (int irange = 0; irange < nrange; ++irange)
    {
        int n_read = INTEGER(VECTOR_ELT(count, 0))[irange];
        SEXP tag = VECTOR_ELT(template_list, TAG_IDX);
        SEXP tmpl;
        if (R_NilValue == tag)
            tmpl = PROTECT(scan_bam_template(R_NilValue));
        else {
            SEXP nms = getAttrib(tag, R_NamesSymbol);
            tmpl = PROTECT(scan_bam_template(nms));
            SEXP count = PROTECT(ScalarInteger(n_read));
            setAttrib(VECTOR_ELT(tmpl, TAG_IDX),
                      install("count"), count);
            UNPROTECT(1);
        }
        for (i = 0; i < LENGTH(names); ++i) {
            if (TAG_IDX == i)
                continue;
            else if (VECTOR_ELT(template_list, i) == R_NilValue)
                SET_VECTOR_ELT(tmpl, i, R_NilValue);
            else {
                if (SEQ_IDX == i) {
                    SET_VECTOR_ELT(tmpl, i, ScalarLogical(TRUE));
                    bd->seq[irange] = new_CharAEAE(n_read, 0);
                } else if (QUAL_IDX == i) {
                    SET_VECTOR_ELT(tmpl, i, ScalarLogical(TRUE));
                    bd->qual[irange] = new_CharAEAE(n_read, 0);
                } else {
                    SEXP elt = allocVector(TYPEOF(VECTOR_ELT(tmpl, i)),
                                           n_read);
                    SET_VECTOR_ELT(tmpl, i, elt);
                }
            }
        }
        SET_VECTOR_ELT(result, irange, tmpl);
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return result;
}

SEXP
_as_XStringSet(CharAEAE *aeae, const char *baseclass)
{
    char classname[40];  /* longest string should be "DNAStringSet" */

    if (snprintf(classname, sizeof(classname), "%sSet", baseclass)
        >= sizeof(classname))
        error("Rsamtools internal error in _as_XStringSet(): "
              "'classname' buffer too small");

    SEXP lkup;
    switch(*baseclass) {
    case 'D':
        lkup = _get_encoding_lookup("B", "DNA");
        break;
    case 'B':
        lkup = R_NilValue;
        break;
    default:
        Rf_error("Rsamtools internal: '%s' unhandled in _as_XStringSet",
                 baseclass);
        break;
    }

    return new_XRawList_from_CharAEAE(classname, baseclass,
                                      aeae, lkup);
}

SEXP
_as_PhredQuality(CharAEAE *aeae)
{
    SEXP xstringset =
        PROTECT(_as_XStringSet(aeae, "BString"));

    SEXP s, t, nmspc, result;
    nmspc = PROTECT(_get_namespace("Rsamtools"));
    NEW_CALL(s, t, "PhredQuality", nmspc, 2);
    CSET_CDR(t, "x", xstringset);
    CEVAL_TO(s, nmspc, result);
    UNPROTECT(2);
    return result;
}

void
_scan_bam_finish1range(_BAM_DATA *bdata, int irange)
{
    SEXP result = VECTOR_ELT((SEXP) bdata->extra, irange), s;
    if ((s = VECTOR_ELT(result, STRAND_IDX)) != R_NilValue) {
        SEXP strand_lvls = PROTECT(_get_strand_levels());
        _as_factor_SEXP(s, strand_lvls);
        UNPROTECT(1);
    }
    bam_header_t *header = bdata->bfile->file->header;
    if ((s = VECTOR_ELT(result, RNAME_IDX)) != R_NilValue) {
        _as_factor(s, (const char **) header->target_name,
                   header->n_targets);
    }
    if ((s = VECTOR_ELT(result, MRNM_IDX)) != R_NilValue) {
        _as_factor(s, (const char **) header->target_name,
                   header->n_targets);
    }
    if ((s = VECTOR_ELT(result, SEQ_IDX)) != R_NilValue)
    {
        PROTECT(s = _as_XStringSet(bdata->seq + irange,
                                   "DNAString"));
        SET_VECTOR_ELT(result, SEQ_IDX, s);
        UNPROTECT(1);
    }
    if ((s = VECTOR_ELT(result, QUAL_IDX)) != R_NilValue)
    {
        PROTECT(s = _as_PhredQuality(bdata->qual + irange));
        SET_VECTOR_ELT(result, QUAL_IDX, s);
        UNPROTECT(1);
    }
}

void
_scan_bam_finish(_BAM_DATA *bdata)
{
    int len = Rf_length((SEXP) bdata->extra);
    for (int irange = 0; irange < len; ++irange)
        _scan_bam_finish1range(bdata, irange);
}

SEXP
_scan_bam(SEXP bfile, SEXP space, SEXP keepFlags, SEXP isSimpleCigar,
	  SEXP filename, SEXP indexname, SEXP filemode,
	  SEXP reverseComplement, SEXP template_list)
{
    SEXP count =
	PROTECT(_count_bam(bfile, space, keepFlags, isSimpleCigar));
    if (R_NilValue == count)
        Rf_error("scanBam failed during countBam");
    if (R_NilValue == space)   /* bam file invalid if fully scanned */
        bamfile_reopen(bfile, filename, indexname, filemode);

    SEXP names = PROTECT(GET_ATTR(template_list, R_NamesSymbol));
    _BAM_DATA *bdata = _init_BAM_DATA(bfile, keepFlags, isSimpleCigar,
				      LOGICAL(reverseComplement)[0]);
    SEXP result =
        PROTECT(_scan_bam_result_init(count, template_list, names,
                                      bdata));
    bdata->extra = (void *) result;

    int status = _do_scan_bam(bdata, space, _scan_bam_parse1);
    if (status < 0) {
        int idx = bdata->idx;
        _BAM_PARSE_STATUS parse_status = bdata->parse_status;
        _Free_BAM_DATA(bdata);
        Rf_error("'scanBam' failed:\n  record: %d\n  error: %d",
                 idx, parse_status);
    }

    if (R_NilValue == space)   /* bam file invalid if fully scanned */
        bamfile_reopen(bfile, filename, indexname, filemode);
    _scan_bam_finish(bdata);
    _Free_BAM_DATA(bdata);
    UNPROTECT(3);
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
    SEXP result = PROTECT(NEW_LIST(2));
    int spc_length =
        (R_NilValue == space) ? 1 : LENGTH(VECTOR_ELT(space, 0));
    SET_VECTOR_ELT(result, 0, NEW_INTEGER(spc_length));
    SET_VECTOR_ELT(result, 1, NEW_NUMERIC(spc_length));
    for (int i = 0; i < spc_length; ++i) {
        INTEGER(VECTOR_ELT(result, 0))[i] =
            REAL(VECTOR_ELT(result, 1))[i] = 0;
    }

    SEXP nms = PROTECT(NEW_CHARACTER(2));
    SET_STRING_ELT(nms, 0, mkChar("records"));
    SET_STRING_ELT(nms, 1, mkChar("nucleotides"));
    setAttrib(result, R_NamesSymbol, nms);
    UNPROTECT(1);

    _BAM_DATA *bdata =
        _init_BAM_DATA(bfile, keepFlags, isSimpleCigar, FALSE);
    bdata->extra = result;

    int status = _do_scan_bam(bdata, space, _count_bam1);
    if (status < 0)
        result = R_NilValue;

    _Free_BAM_DATA(bdata);
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
    _BAM_DATA *bdata =
        _init_BAM_DATA(bfile, keepFlags, isSimpleCigar, FALSE);
    /* FIXME: this just copies the header... */
    bam_header_t *header = BAMFILE(bfile)->file->header;
    samfile_t *f_out =
        _bam_tryopen(translateChar(STRING_ELT(fout_name, 0)),
                     CHAR(STRING_ELT(fout_mode, 0)), header);
    bdata->extra = f_out;

    int status = _do_scan_bam(bdata, space, _filter_bam1);

    /* sort and index destintation ? */
    /* cleanup */
    samclose(f_out);
    _Free_BAM_DATA(bdata);

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
