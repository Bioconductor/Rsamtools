#include "samtools/khash.h"
#include "samtools/sam.h"
#include "bam_data.h"
#include "scan_bam_data.h"
#include "io_sam.h"
#include "bamfile.h"
#include "encode.h"
#include "utilities.h"
#include "S4Vectors_interface.h"
#include "XVector_interface.h"
#include "Biostrings_interface.h"
#include "bam_mate_iter.h"

/* from samtoools/bam_sort.c */
void bam_sort_core(int is_by_qname, const char *fn, const char *prefix,
                   size_t max_mem);

#define SEQUENCE_BUFFER_ALLOCATION_ERROR 1

static const char *TMPL_ELT_NMS[] = {
    "qname", "flag", "rname", "strand", "pos", "qwidth", "mapq", "cigar",
    "mrnm", "mpos", "isize", "seq", "qual", "tag", "groupid", "mate_status"
    /* "vtype", "value" */
};

static const int N_TMPL_ELTS = sizeof(TMPL_ELT_NMS) / sizeof(const char *);

/* utility */

void _check_is_bam(const char *filename)
{
    int magic_len;
    char buf[4];
    bamFile bfile = bam_open(filename, "r");

    if (bfile == 0)
        Rf_error("failed to open SAM/BAM file\n  file: '%s'", filename);

    magic_len = bam_read(bfile, buf, 4);
    bam_close(bfile);

    if (magic_len != 4 || strncmp(buf, "BAM\001", 4) != 0)
        Rf_error("'filename' is not a BAM file\n  file: %s", filename);
}

/* template */

void _bam_check_template_list(SEXP template_list)
{
    if (!IS_LIST(template_list) || LENGTH(template_list) != N_TMPL_ELTS)
        Rf_error("'template' must be list(%d)", N_TMPL_ELTS);
    SEXP names = GET_ATTR(template_list, R_NamesSymbol);
    if (!IS_CHARACTER(names) || LENGTH(names) != N_TMPL_ELTS)
        Rf_error("'names(template)' must be character(%d)", N_TMPL_ELTS);
    for (int i = 0; i < LENGTH(names); ++i)
        if (strcmp(TMPL_ELT_NMS[i], CHAR(STRING_ELT(names, i))) != 0)
            Rf_error("'template' names do not match scan_bam_template\n'");
}

static SEXP _tmpl_strand()
{
    SEXP strand = PROTECT(NEW_INTEGER(0));
    _as_strand(strand);
    UNPROTECT(1);
    return strand;
}

static SEXP _tmpl_DNAStringSet()
{
    CharAEAE *aeae = new_CharAEAE(0, 0);
    SEXP lkup = PROTECT(_get_lkup("DNAString"));
    SEXP ans = new_XRawList_from_CharAEAE("DNAStringSet", "DNAString",
                                          aeae, lkup);
    UNPROTECT(1);
    return ans;
}

static SEXP _tmpl_BStringSet()
{
    CharAEAE *aeae = new_CharAEAE(0, 0);
    return new_XRawList_from_CharAEAE("BStringSet", "BString", aeae,
                                      R_NilValue);
}

static SEXP _tmpl_PhredQuality()
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

SEXP scan_bam_template(SEXP rname, SEXP tag)
{
    if (R_NilValue != tag)
        if (!IS_CHARACTER(tag))
            Rf_error("'tag' must be NULL or 'character()'");
    SEXP tmpl = PROTECT(NEW_LIST(N_TMPL_ELTS));
    SET_VECTOR_ELT(tmpl, QNAME_IDX, NEW_CHARACTER(0));
    SET_VECTOR_ELT(tmpl, FLAG_IDX, NEW_INTEGER(0));
    SET_VECTOR_ELT(tmpl, RNAME_IDX, rname);
    SET_VECTOR_ELT(tmpl, STRAND_IDX, _tmpl_strand());
    SET_VECTOR_ELT(tmpl, POS_IDX, NEW_INTEGER(0));
    SET_VECTOR_ELT(tmpl, QWIDTH_IDX, NEW_INTEGER(0));
    SET_VECTOR_ELT(tmpl, MAPQ_IDX, NEW_INTEGER(0));
    SET_VECTOR_ELT(tmpl, CIGAR_IDX, NEW_CHARACTER(0));
    SET_VECTOR_ELT(tmpl, MRNM_IDX, rname);
    SET_VECTOR_ELT(tmpl, MPOS_IDX, NEW_INTEGER(0));
    SET_VECTOR_ELT(tmpl, ISIZE_IDX, NEW_INTEGER(0));
    SET_VECTOR_ELT(tmpl, SEQ_IDX, _tmpl_DNAStringSet());
    SET_VECTOR_ELT(tmpl, QUAL_IDX, _tmpl_PhredQuality());
    SET_VECTOR_ELT(tmpl, PARTITION_IDX, NEW_INTEGER(0));
    SET_VECTOR_ELT(tmpl, MATES_IDX, NEW_INTEGER(0));
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

/* header */

SEXP _read_bam_header(SEXP ext, SEXP what)
{
    samfile_t *sfile = BAMFILE(ext)->file;
    bam_header_t *header = sfile->header;

    SEXP ans = PROTECT(NEW_LIST(2));
    SEXP nms = NEW_CHARACTER(2);
    setAttrib(ans, R_NamesSymbol, nms);
    SET_STRING_ELT(nms, 0, mkChar("targets"));
    SET_STRING_ELT(nms, 1, mkChar("text"));

    if (LOGICAL(what)[0] == TRUE) { /* 'targets' */
        int n_elts = header->n_targets;
        SET_VECTOR_ELT(ans, 0, NEW_INTEGER(n_elts));
        SEXP tlen = VECTOR_ELT(ans, 0);   /* target length */
        SEXP tnm = NEW_CHARACTER(n_elts); /* target name */
        setAttrib(tlen, R_NamesSymbol, tnm);
        for (int j = 0; j < n_elts; ++j) {
            INTEGER(tlen)[j] = header->target_len[j];
            SET_STRING_ELT(tnm, j, mkChar(header->target_name[j]));
        }
    }

    if (LOGICAL(what)[1] == TRUE) { /* 'text' */
        int n_text_elts = 0;
        for (int i = 0; i < header->l_text; ++i)
            if (header->text[i] == '\n')
                n_text_elts += 1;
        SET_VECTOR_ELT(ans, 1, NEW_LIST(n_text_elts));
        SEXP text = VECTOR_ELT(ans, 1);
        SEXP tag = NEW_CHARACTER(n_text_elts);
        setAttrib(text, R_NamesSymbol, tag);

        int start = 0, end;
        for (int i = 0; i < n_text_elts; ++i) {
            int n_elts = header->text[start] == '\n' ? 0 : 1;
            end = start;
            while (header->text[end] != '\n') {
                if (header->text[end] == '\t')
                    ++n_elts;
                ++end;
            }
            if (n_elts == 0) {
                SET_VECTOR_ELT(text, i, NEW_CHARACTER(0));
                /* SET_STRING_ELT(tag, i, mkChar("")); */
                start = end + 1;
                continue;
            }
            SET_VECTOR_ELT(text, i, NEW_CHARACTER(n_elts - 1));
            SEXP elts = VECTOR_ELT(text, i);

            for (int j = 0; j < n_elts; ++j) {
                end = start;
                while (header->text[end] != '\t' && header->text[end] != '\n')
                    ++end;
                SEXP elt = mkCharLen(&header->text[start], end - start);
                if (j == 0)   /* tag */
                    SET_STRING_ELT(tag, i, elt);
                else
                    SET_STRING_ELT(elts, j - 1, elt);
                start = end + 1;
            }
        }
    }

    UNPROTECT(1);
    return ans;
}

/* scan_bam framework */

int check_qname(char *last_qname, int bufsize, bam1_t *bam, int max) 
{
    if (0 != strcmp(last_qname, bam1_qname(bam))) {
        /* stop reading */
        if (max) {
            return -1;
        /* get next qname, continue reading */
        } else {
            if (bam->core.l_qname > bufsize) {
                Free(last_qname);
                bufsize = bam->core.l_qname;
                last_qname = Calloc(bufsize, char);
            }
            strcpy(last_qname, bam1_qname(bam));
            return 1;
        }
    /* same qname, continue reading */
    } else {
        return 0;
    }
}

int _samread(BAM_FILE bfile, BAM_DATA bd, const int yieldSize,
             bam_fetch_f parse1)
{
    int yield = 0, status = 1, bufsize = 1000;
    char *last_qname = Calloc(bufsize, char);
    bam1_t *bam = bam_init1();

    while (samread(bfile->file, bam) >= 0) {
        if (NA_INTEGER != yieldSize) {
            if (bd->obeyQname)
                status = check_qname(last_qname, bufsize, bam, 
                                     yield >= yieldSize);
            if (status < 0)
                break;
        }
 
        int result = parse1(bam, bd);
        if (result < 0) {   /* parse error: e.g., cigar buffer overflow */
            bam_destroy1(bam);
            Free(last_qname);
            return yield;
        } else if (result == 0L) /* does not pass filter */
            continue;

        yield += status;
        if (NA_INTEGER != yieldSize && yield == yieldSize) { 
            bfile->pos0 = bam_tell(bfile->file->x.bam);
            if (!bd->obeyQname) 
                break;
        }
    }

    bam_destroy1(bam);
    Free(last_qname);
    return yield;
}

int _samread_mate(BAM_FILE bfile, BAM_DATA bd, const int yieldSize,
                  bam_fetch_mate_f parse1_mate)
{
    int yield = 0;
    bam_mates_t *bam_mates = bam_mates_new();

    while (samread_mate(bfile->file->x.bam, bfile->index,
                        &bfile->iter, bam_mates, bd) > 0) {

        if (NA_INTEGER != yieldSize && yield  >= yieldSize)
            break;

        int result = parse1_mate(bam_mates, bd);
        if (result < 0) {       /* parse error */
            bam_mates_destroy(bam_mates);
            return result;
        } else if (result == 0)
            continue;

        yield += 1;
        if (NA_INTEGER != yieldSize && yield == yieldSize) { 
            bfile->pos0 = bam_tell(bfile->file->x.bam);
            break;
        }

    }

    bam_mates_destroy(bam_mates);
    return yield;
}

/* read complete file */
static int _scan_bam_all(BAM_DATA bd, bam_fetch_f parse1,
                         bam_fetch_mate_f parse1_mate, _FINISH1_FUNC finish1)
{
    BAM_FILE bfile = _bam_file_BAM_DATA(bd);
    const int yieldSize = bd->yieldSize;
    int yield = 0;

    bam_seek(bfile->file->x.bam, bfile->pos0, SEEK_SET);
    if (bd->asMates) {
        yield = _samread_mate(bfile, bd, yieldSize, parse1_mate);
    } else {
        yield = _samread(bfile, bd, yieldSize, parse1);
    }

    /* end-of-file */
    if ((NA_INTEGER == yieldSize) || (yield < yieldSize))
        bfile->pos0 = bam_tell(bfile->file->x.bam);
    if ((NULL != finish1) && (bd->iparsed >= 0))
        (*finish1) (bd);

    return bd->iparsed;
}

/* read ranges */
static int _scan_bam_fetch(BAM_DATA bd, SEXP space, int *start, int *end,
                           bam_fetch_f parse1, bam_fetch_mate_f parse1_mate,
                           _FINISH1_FUNC finish1)
{
    int tid;
    BAM_FILE bfile = _bam_file_BAM_DATA(bd);
    samfile_t *sfile = bfile->file;
    bam_index_t *bindex = bfile->index;
    const int initial = bd->iparsed;

    for (int irange = bfile->irange0; irange < LENGTH(space); ++irange) {
        const char *spc = translateChar(STRING_ELT(space, irange));
        const int starti =
            start[irange] > 0 ? start[irange] - 1 : start[irange];
        for (tid = 0; tid < sfile->header->n_targets; ++tid) {
            if (strcmp(spc, sfile->header->target_name[tid]) == 0)
                break;
        }
        if (tid == sfile->header->n_targets) {
            Rf_warning("space '%s' not in BAM header", spc);
            bd->irange += 1;
            return -1;
        }
        if (bd->asMates) {
            bam_fetch_mate(sfile->x.bam, bindex, tid, starti, end[irange], 
                           bd, parse1_mate);
        } else {
            bam_fetch(sfile->x.bam, bindex, tid, starti, end[irange],
                      bd, parse1);
        }

        if (NULL != finish1)
            (*finish1) (bd);
        bd->irange += 1;
        if ((NA_INTEGER != bd->yieldSize) &&
            (bd->iparsed - initial >= bd->yieldSize))
            break;
    }
    bfile->irange0 = bd->irange;
        
    return bd->iparsed - initial;
}

int _do_scan_bam(BAM_DATA bd, SEXP space, bam_fetch_f parse1,
                 bam_fetch_mate_f parse1_mate, _FINISH1_FUNC finish1)
{
    int status;

    if (R_NilValue == space)
        /* everything */
        status = _scan_bam_all(bd, parse1, parse1_mate, finish1);
    else {                      
        /* fetch */
        BAM_FILE bfile = _bam_file_BAM_DATA(bd);
        if (NULL == bfile->index)
            Rf_error("valid 'index' file required");
        status = _scan_bam_fetch(bd, VECTOR_ELT(space, 0),
                                 INTEGER(VECTOR_ELT(space, 1)),
                                 INTEGER(VECTOR_ELT(space, 2)),
                                 parse1, parse1_mate, finish1);
    }

    return status;
}

/* scan_bam */

static int _filter_and_parse1(const bam1_t *bam, void *data)
{
    /* requires 'data' to be BAM_DATA */
    BAM_DATA bd = (BAM_DATA) data; 
    int result = _filter_and_parse1_BAM_DATA(bam, bd);
    if (result < 0) {
        /* parse error: e.g., cigar buf overflow */
        _grow_SCAN_BAM_DATA(bd, 0);
        bd->iparsed = -1;
    }
    return result;
}

static int _filter_and_parse1_mate(const bam_mates_t *mates, void *data)
{
    /* requires 'data' to be BAM_DATA, data->extra to be SCAN_BAM_DATA */
    BAM_DATA bd = (BAM_DATA) data;
    SCAN_BAM_DATA sbd = (SCAN_BAM_DATA) bd->extra;
    int yield = 0, pass = 0;

    sbd->mates_flag = mates->mated == MATE_UNKNOWN ? NA_INTEGER : mates->mated;
    sbd->partition_id += 1;

    for (int i = 0; i < mates->n; ++i) {
        int result = _filter_and_parse1_BAM_DATA(mates->bams[i], bd);
        if (result < 0) {
            _grow_SCAN_BAM_DATA(bd, 0);
            return result;
        }
        pass += result;
    }

    if (pass > 0) {
        yield = 1;
    } else {
        sbd->partition_id -= 1;
    }

    return yield;
}

SEXP _scan_bam_result_init(SEXP template_list, SEXP names, SEXP space,
                           BAM_FILE bfile)
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
    bam_header_t *header = bfile->file->header;
    SEXP rname = PROTECT(NEW_INTEGER(0));
    _as_factor(rname, (const char **) header->target_name, header->n_targets);
    
    for (int irange = 0; irange < nrange; ++irange) {
        SEXP tag = VECTOR_ELT(template_list, TAG_IDX);
        SEXP tmpl;
        if (R_NilValue == tag)
            tmpl = PROTECT(scan_bam_template(rname, R_NilValue));
        else {
            SEXP nms = getAttrib(tag, R_NamesSymbol);
            tmpl = PROTECT(scan_bam_template(rname, nms));
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
    UNPROTECT(2);
    return result;
}

SEXP _scan_bam(SEXP bfile, SEXP space, SEXP keepFlags, SEXP isSimpleCigar,
               SEXP tagFilter, SEXP reverseComplement, SEXP yieldSize,
               SEXP template_list, SEXP obeyQname, SEXP asMates,
               SEXP qnamePrefixEnd, SEXP qnameSuffixStart)
{
    SEXP names = PROTECT(GET_ATTR(template_list, R_NamesSymbol));
    SEXP result = PROTECT(_scan_bam_result_init(template_list, names, space,
                                                BAMFILE(bfile)));
    SCAN_BAM_DATA sbd = _init_SCAN_BAM_DATA(result);

    char qname_prefix = '\0';
    SEXP prefix_elt = STRING_ELT(qnamePrefixEnd, 0);
    if (prefix_elt != NA_STRING)
        qname_prefix = CHAR(prefix_elt)[0];
    char qname_suffix = '\0';
    SEXP suffix_elt = STRING_ELT(qnameSuffixStart, 0);
    if (suffix_elt != NA_STRING)
        qname_suffix = CHAR(suffix_elt)[0];

    BAM_DATA bd = _init_BAM_DATA(bfile, space, keepFlags, isSimpleCigar,
                                 tagFilter,
                                 LOGICAL(reverseComplement)[0],
                                 INTEGER(yieldSize)[0],
                                 LOGICAL(obeyQname)[0], 
                                 LOGICAL(asMates)[0], 
                                 qname_prefix, qname_suffix, (void *) sbd);

    int status = _do_scan_bam(bd, space, _filter_and_parse1,
                              _filter_and_parse1_mate, _finish1range_BAM_DATA);
    if (status < 0) {
        int idx = bd->irec;
        int parse_status = bd->parse_status;
        _Free_SCAN_BAM_DATA(bd->extra);
        _Free_BAM_DATA(bd);
        Rf_error("'scanBam' failed:\n  record: %d\n  error: %d",
                 idx, parse_status);
    }

    _Free_SCAN_BAM_DATA(bd->extra);
    _Free_BAM_DATA(bd);
    UNPROTECT(2);
    return result;
}

/* count */

static int _count1(const bam1_t * bam, void *data)
{
    return _count1_BAM_DATA(bam, (BAM_DATA) data);
}

SEXP _count_bam(SEXP bfile, SEXP space, SEXP keepFlags, SEXP isSimpleCigar,
    SEXP tagFilter)
{
    SEXP result = PROTECT(NEW_LIST(2));
    BAM_DATA bd =
        _init_BAM_DATA(bfile, space, keepFlags, isSimpleCigar, tagFilter, 0,
                       NA_INTEGER, 0, 0, '\0', '\0', result);

    SET_VECTOR_ELT(result, 0, NEW_INTEGER(bd->nrange));
    SET_VECTOR_ELT(result, 1, NEW_NUMERIC(bd->nrange));
    for (int i = 0; i < bd->nrange; ++i) {
        INTEGER(VECTOR_ELT(result, 0))[i] = REAL(VECTOR_ELT(result, 1))[i] = 0;
    }

    SEXP nms = PROTECT(NEW_CHARACTER(2));
    SET_STRING_ELT(nms, 0, mkChar("records"));
    SET_STRING_ELT(nms, 1, mkChar("nucleotides"));
    setAttrib(result, R_NamesSymbol, nms);
    UNPROTECT(1);

    int status = _do_scan_bam(bd, space, _count1, NULL, NULL);
    if (status < 0) {
        int idx = bd->irec;
        int parse_status = bd->parse_status;
        _Free_BAM_DATA(bd);
        UNPROTECT(1);
        Rf_error("'countBam' failed:\n  record: %d\n  error: %d",
                 idx, parse_status);
    }

    _Free_BAM_DATA(bd);
    UNPROTECT(1);
    return result;
}

void scan_bam_cleanup()
{
    /* placeholder */
}

/* filterBam */

static int _prefilter1(const bam1_t *bam, void *data)
{
    BAM_DATA bd = (BAM_DATA) data;
    bd->irec += 1;
    if (!_filter1_BAM_DATA(bam, bd))
        return 0;
    bambuffer_push((BAM_BUFFER) bd->extra, bam);
    bd->iparsed += 1;
    return 1;
}

static int _prefilter1_mate(const bam_mates_t *mates, void *data)
{
    BAM_DATA bd = (BAM_DATA) data;
    BAM_BUFFER buf = (BAM_BUFFER) bd->extra;
    int parsed;

    buf->partition_id += 1;
    buf->mate_flag = mates->mated == MATE_UNKNOWN ? NA_INTEGER : mates->mated;
    parsed = 0;
    for (int i = 0; i < mates->n; ++i)
        parsed += _prefilter1(mates->bams[i], data);

    if (parsed == 0)
        buf->partition_id -= 1;

    return parsed;
}

SEXP
_prefilter_bam(SEXP bfile, SEXP space, SEXP keepFlags, SEXP isSimpleCigar,
               SEXP tagFilter, SEXP yieldSize, SEXP obeyQname, SEXP asMates, 
               SEXP qnamePrefixEnd, SEXP qnameSuffixStart)
{
    SEXP ext = PROTECT(bambuffer(INTEGER(yieldSize)[0],
                                 LOGICAL(asMates)[0]));
    char qname_prefix = '\0';
    SEXP prefix_elt = STRING_ELT(qnamePrefixEnd, 0);
    if (prefix_elt != NA_STRING)
        qname_prefix = CHAR(prefix_elt)[0];
    char qname_suffix = '\0';
    SEXP suffix_elt = STRING_ELT(qnameSuffixStart, 0);
    if (suffix_elt != NA_STRING)
        qname_suffix = CHAR(suffix_elt)[0];
    BAM_DATA bd = _init_BAM_DATA(bfile, space, keepFlags, isSimpleCigar,
                                 tagFilter, 0, INTEGER(yieldSize)[0],
                                 LOGICAL(obeyQname)[0], 
                                 LOGICAL(asMates)[0], 
                                 qname_prefix, qname_suffix, BAMBUFFER(ext));
    int status =
        _do_scan_bam(bd, space, _prefilter1, _prefilter1_mate, NULL);
    if (status < 0) {
        int idx = bd->irec;
        int parse_status = bd->parse_status;
        _Free_BAM_DATA(bd);
        UNPROTECT(1);
        Rf_error("'filterBam' prefilter failed:\n  record: %d\n  error: %d",
                 idx, parse_status);
    }

    _Free_BAM_DATA(bd);
    UNPROTECT(1);
    return ext;
}

static int _filter1(const bam1_t * bam, void *data)
{
    BAM_DATA bd = (BAM_DATA) data;
    bd->irec += 1;
    if (!_filter1_BAM_DATA(bam, bd))
        return 0;
    samwrite((samfile_t *) bd->extra, bam);
    bd->iparsed += 1;
    return 1;
}

SEXP
_filter_bam(SEXP bfile, SEXP space, SEXP keepFlags,
            SEXP isSimpleCigar, SEXP tagFilter, SEXP fout_name, SEXP fout_mode)
{
    /* open destination */
    BAM_DATA bd =
        _init_BAM_DATA(bfile, space, keepFlags, isSimpleCigar, tagFilter, 0,
                       NA_INTEGER, 0, 0, '\0', '\0', NULL);
    /* FIXME: this just copies the header... */
    bam_header_t *header = BAMFILE(bfile)->file->header;
    samfile_t *f_out = _bam_tryopen(translateChar(STRING_ELT(fout_name, 0)),
                                    CHAR(STRING_ELT(fout_mode, 0)), header);
    bd->extra = f_out;

    int status = _do_scan_bam(bd, space, _filter1, NULL, NULL);
    if (status < 0) {
        int idx = bd->irec;
        int parse_status = bd->parse_status;
        _Free_BAM_DATA(bd);
        samclose(f_out);
        Rf_error("'filterBam' failed:\n  record: %d\n  error: %d",
                 idx, parse_status);
    }

    /* sort and index destintation ? */
    /* cleanup */
    _Free_BAM_DATA(bd);
    samclose(f_out);

    return status < 0 ? R_NilValue : fout_name;
}

/* merge_bam */

/* from bam_sort.c */

#define MERGE_RG     1
#define MERGE_LEVEL1 4
#define MERGE_FORCE  8

int bam_merge_core(int by_qname, const char *out, const char *headers,
                   int n, char * const *fn, int flag, const char *reg);

SEXP merge_bam(SEXP fnames, SEXP destination, SEXP overwrite,
               SEXP hname, SEXP regionStr, SEXP isByQname,
               SEXP addRG, SEXP compressLevel1)
{
    int i;

    if (!IS_CHARACTER(fnames) || 2 > Rf_length(fnames))
        Rf_error("'files' must be a character() with length >= 2");
    if (!IS_CHARACTER(hname) || 1 <  Rf_length(hname))
        Rf_error("'header' must be character() with length <= 1");
    if (!IS_CHARACTER(destination) || 1 != Rf_length(destination))
        Rf_error("'destination' must be character(1)");
    if (!IS_LOGICAL(overwrite) || 1 != Rf_length(overwrite))
        Rf_error("'overwrite' must be logical(1)");
    if (!IS_CHARACTER(regionStr) || 1 <  Rf_length(regionStr))
        Rf_error("'region' must define 0 or 1 regions");
    if (!IS_LOGICAL(isByQname) || 1 != Rf_length(isByQname))
        Rf_error("'isByQname' must be logical(1)");
    if (!IS_LOGICAL(addRG) || 1 != Rf_length(addRG))
        Rf_error("'addRG' must be logical(1)");
    if (!IS_LOGICAL(compressLevel1) || 1 != Rf_length(compressLevel1))
        Rf_error("'compressLevel1' must be logical(1)");

    char ** fileNames = (char **)
        R_alloc(sizeof(const char *), Rf_length(fnames));
    for (i = 0; i < Rf_length(fnames); ++i)
        fileNames[i] = (char *) translateChar(STRING_ELT(fnames, i));

    const char *hfName = 0 == Rf_length(hname) ?
        NULL : translateChar(STRING_ELT(hname, 0));

    int flag = 0;
    if (LOGICAL(addRG)[0])
        flag = ((int) flag) | ((int) MERGE_RG);
    if (LOGICAL(overwrite)[0])
        flag = ((int) flag) | ((int) MERGE_FORCE);
    if (LOGICAL(compressLevel1)[0])
        flag = ((int) flag) | ((int) MERGE_LEVEL1);

    const char *region = 0 == Rf_length(regionStr) ?
        NULL : translateChar(STRING_ELT(regionStr, 0));

    int res = bam_merge_core(LOGICAL(isByQname)[0],
                             translateChar(STRING_ELT(destination, 0)),
                             hfName, Rf_length(fnames), fileNames,
                             flag, region);
    if (res < 0)
        Rf_error("'mergeBam' failed with error code %d", res);

    return destination;
}

/* sort_bam */

SEXP sort_bam(SEXP filename, SEXP destination, SEXP isByQname, SEXP maxMemory)
{
    if (!IS_CHARACTER(filename) || 1 != LENGTH(filename))
        Rf_error("'filename' must be character(1)");
    if (!IS_CHARACTER(destination) || 1 != LENGTH(destination))
        Rf_error("'destination' must be character(1)");
    if (!IS_LOGICAL(isByQname) || LENGTH(isByQname) != 1)
        Rf_error("'isByQname' must be logical(1)");
    if (!IS_INTEGER(maxMemory) || LENGTH(maxMemory) != 1 ||
        INTEGER(maxMemory)[0] < 1)
        Rf_error("'maxMemory' must be a positive integer(1)");

    const char *fbam = translateChar(STRING_ELT(filename, 0));
    const char *fout = translateChar(STRING_ELT(destination, 0));
    int sortMode = asInteger(isByQname);

    size_t maxMem = (size_t) INTEGER(maxMemory)[0] * 1024 * 1024;
    _check_is_bam(fbam);
    bam_sort_core(sortMode, fbam, fout, maxMem);

    return destination;
}

/* index_bam */

SEXP index_bam(SEXP indexname)
{
    if (!IS_CHARACTER(indexname) || 1 != LENGTH(indexname))
        Rf_error("'indexname' must be character(1)");
    const char *fbam = translateChar(STRING_ELT(indexname, 0));

    _check_is_bam(fbam);
    int status = bam_index_build(fbam);

    if (0 != status)
        Rf_error("failed to build index\n  file: %s", fbam);
    char *fidx = (char *) R_alloc(strlen(fbam) + 5, sizeof(char));
    sprintf(fidx, "%s.bai", fbam);
    return mkString(fidx);
}
