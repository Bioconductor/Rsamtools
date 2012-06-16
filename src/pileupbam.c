#include "samtools/bam.h"
#include "samtools/khash.h"
#include "pileupbam.h"
#include "bamfile.h"
#include "utilities.h"

typedef enum {
    YIELDBY_RANGE = 0, YIELDBY_POSITION
} YIELDBY;

typedef enum {
    WHAT_SEQ = 1, WHAT_QUAL = 2
} WHAT;

typedef struct {
    _BAM_FILE *bfile;
    bamFile fp;
    bam_iter_t iter;
    /* read filter params */
    int min_map_quality;
    uint32_t keep_flag[2];
} BAM_ITER_T;

typedef struct {
    int n_files, *n_plp;
    BAM_ITER_T **mfile;
    const bam_pileup1_t **plp;
    bam_mplp_t mplp_iter;
} PILEUP_ITER_T;

typedef struct {
    const char *chr;
    int i_spc, start, end;
} SPACE_T;

typedef struct {
    SPACE_T *spc;
    int i_spc, n_spc;
    const char **chr;
    int *start, *end, stashed;
} SPACE_ITER_T;

typedef struct {
    int i_yld;
    int *pos, *seq, *qual;
} PILEUP_RESULT_T;

typedef struct {
    int n_files;
    SEXP names;
    int min_base_quality, min_map_quality, min_depth, max_depth;
    uint32_t keep_flag[2];
    int yieldSize, yieldAll;
    YIELDBY yieldBy;
    WHAT what;
} PILEUP_PARAM_T;

/* from bam_aux.c; should really be exported part of library? */
KHASH_MAP_INIT_STR(s, int)

const int QUAL_LEVELS = 94,     /* printable ASCII */
    SEQ_LEVELS = 5;             /* A, C, G, T, N */

static void _bam_header_hash_init(bam_header_t * header)
{
    if (0 == header->hash) {
        int ret, i;
        khiter_t iter;
        khash_t(s) * h;
        header->hash = h = kh_init(s);
        for (i = 0; i < header->n_targets; ++i) {
            iter = kh_put(s, h, header->target_name[i], &ret);
            kh_value(h, iter) = i;
        }
    }
}

/* PILEUP_ITER_T */

static PILEUP_ITER_T *_iter_init(SEXP files, PILEUP_PARAM_T * param)
{
    int i;
    PILEUP_ITER_T *iter = Calloc(1, PILEUP_ITER_T);
    iter->n_files = Rf_length(files);;
    iter->mfile = Calloc(iter->n_files, BAM_ITER_T *);
    iter->mfile[0] = Calloc(iter->n_files, BAM_ITER_T);
    for (i = 0; i < iter->n_files; ++i) {
        iter->mfile[i] = iter->mfile[0] + i;
        iter->mfile[i]->bfile = BAMFILE(VECTOR_ELT(files, i));
        iter->mfile[i]->fp = iter->mfile[i]->bfile->file->x.bam;
        iter->mfile[i]->min_map_quality = param->min_map_quality;
        iter->mfile[i]->keep_flag[0] = param->keep_flag[0];
        iter->mfile[i]->keep_flag[1] = param->keep_flag[1];
        /* header hash destroyed when file closed */
        _bam_header_hash_init(iter->mfile[i]->bfile->file->header);
    }

    iter->plp = Calloc(iter->n_files, const bam_pileup1_t *);
    iter->n_plp = Calloc(iter->n_files, int);

    return iter;
}

static PILEUP_ITER_T *_iter_destroy(PILEUP_ITER_T * iter)
{
    Free(iter->plp);
    Free(iter->n_plp);
    Free(iter->mfile[0]);
    Free(iter->mfile);
    Free(iter);
    return NULL;
}

/* SPACE_ITER_T */

static SPACE_ITER_T *_space_iter_init(SEXP space)
{
    int i;
    SPACE_ITER_T *iter;

    iter = Calloc(1, SPACE_ITER_T);

    iter->i_spc = -1;
    iter->n_spc = Rf_length(VECTOR_ELT(space, 0));
    iter->chr = Calloc(iter->n_spc, const char *);
    for (i = 0; i < iter->n_spc; ++i)
        iter->chr[i] = CHAR(STRING_ELT(VECTOR_ELT(space, 0), i));
    iter->start = INTEGER(VECTOR_ELT(space, 1));
    iter->end = INTEGER(VECTOR_ELT(space, 2));

    iter->stashed = FALSE;
    iter->spc = Calloc(1, SPACE_T);

    return iter;
}

static SPACE_T *_space_iter_next(SPACE_ITER_T * iter)
{
    iter->stashed = FALSE;
    iter->i_spc += 1;
    SPACE_T *spc = iter->spc;
    if (iter->i_spc < iter->n_spc) {
        spc->i_spc = iter->i_spc;
        spc->chr = iter->chr[iter->i_spc];
        spc->start = iter->start[iter->i_spc];
        spc->end = iter->end[iter->i_spc];
    } else {
        spc = NULL;
    }
    return spc;
}

static void _space_iter_stash(SPACE_ITER_T * iter, SPACE_T * spc)
{
    iter->spc->chr = spc->chr;
    iter->spc->start = spc->start;
    iter->spc->end = spc->end;
    iter->stashed = TRUE;
}

static SPACE_T *_space_iter_pop(SPACE_ITER_T * iter)
{
    SPACE_T *spc = NULL;
    if (iter->stashed) {
        spc = iter->spc;
        iter->stashed = FALSE;
    }
    return spc;
}

static SPACE_ITER_T *_space_iter_destroy(SPACE_ITER_T * iter)
{
    Free(iter->chr);
    Free(iter->spc);
    Free(iter);
    return NULL;
}

/*  */

static SEXP _lst_elt(SEXP lst, const char *name, const char *lst_name)
{
    SEXP nms = GET_NAMES(lst);
    SEXP elt_nm = PROTECT(mkChar(name));
    int i;
    for (i = 0; i < Rf_length(nms); ++i)
        if (elt_nm == STRING_ELT(nms, i))
            break;
    UNPROTECT(1);
    if (i == Rf_length(nms))
        Rf_error("'%s' does not contain element '%s'", lst_name, name);
    return VECTOR_ELT(lst, i);
}

static int _mplp_read_bam(void *data, bam1_t * b)
{
    BAM_ITER_T *mdata = (BAM_ITER_T *) data;
    uint32_t test_flag;
    int skip, result;

    do {
        result = mdata->iter ?
            bam_iter_read(mdata->fp, mdata->iter, b) : bam_read1(mdata->fp, b);
        if (0 >= result)
            break;

        skip = FALSE;
        test_flag = (mdata->keep_flag[0] & ~b->core.flag) |
            (mdata->keep_flag[1] & b->core.flag);
        if (~test_flag & 2047u)
            skip = TRUE;
        else if (b->core.tid < 0 || (b->core.flag & BAM_FUNMAP))
            skip = TRUE;
        else if (b->core.qual < mdata->min_map_quality)
            skip = TRUE;
    } while (skip);
    return result;
}

static SEXP _seq_rle(int *cnt, const char **chr, int n)
{
    int i = 0, j;
    SEXP s, t;

    for (j = 1; j < n; ++j) {
        if (0 != strcmp(chr[j], chr[j - 1])) {
            i += 1;
            chr[i] = chr[j];
            cnt[i] = cnt[j] - cnt[i - 1];
        } else
            cnt[i] += cnt[j] - cnt[j - 1];
    }
    if (n)
        n = i + 1;

    s = PROTECT(NEW_INTEGER(n));
    t = NEW_CHARACTER(n);
    Rf_setAttrib(s, R_NamesSymbol, t);

    for (i = 0; i < n; ++i) {
        INTEGER(s)[i] = cnt[i];
        SET_STRING_ELT(t, i, mkChar(chr[i]));
    }

    UNPROTECT(1);
    return s;
}

static SEXP _mplp_setup_R(const PILEUP_PARAM_T * param,
                          PILEUP_RESULT_T * result)
{
#ifdef PILEUPBAM_DEBUG
    REprintf("_mplp_setup_R\n");
#endif
    SEXP alloc = PROTECT(NEW_LIST(4)),
        nms = PROTECT(NEW_CHARACTER(4)), opos, oseq, oqual, dimnms, dim_elt;
    char qualbuf[] = { ' ', '\0' };
    int i;

    SET_STRING_ELT(nms, 0, mkChar("seqnames"));
    SET_STRING_ELT(nms, 1, mkChar("pos"));
    SET_STRING_ELT(nms, 2, mkChar("seq"));
    SET_STRING_ELT(nms, 3, mkChar("qual"));
    Rf_setAttrib(alloc, R_NamesSymbol, nms);

    result->i_yld = 0;

    SET_VECTOR_ELT(alloc, 0, _seq_rle(NULL, NULL, 0));

    opos = NEW_INTEGER(param->yieldSize);
    memset(INTEGER(opos), 0, sizeof(int) * Rf_length(opos));
    SET_VECTOR_ELT(alloc, 1, opos);
    result->pos = INTEGER(opos);

    if (param->what & WHAT_SEQ) {
        oseq = Rf_alloc3DArray(INTSXP, SEQ_LEVELS, param->n_files,
                               param->yieldSize);
        memset(INTEGER(oseq), 0, sizeof(int) * Rf_length(oseq));
        SET_VECTOR_ELT(alloc, 2, oseq);

        dimnms = NEW_LIST(3);
        Rf_setAttrib(oseq, R_DimNamesSymbol, dimnms);

        dim_elt = NEW_CHARACTER(SEQ_LEVELS);
        SET_VECTOR_ELT(dimnms, 0, dim_elt);
        SET_VECTOR_ELT(dimnms, 1, param->names);
        SET_VECTOR_ELT(dimnms, 2, R_NilValue);

        SET_STRING_ELT(dim_elt, 0, mkChar("A"));
        SET_STRING_ELT(dim_elt, 1, mkChar("C"));
        SET_STRING_ELT(dim_elt, 2, mkChar("G"));
        SET_STRING_ELT(dim_elt, 3, mkChar("T"));
        SET_STRING_ELT(dim_elt, 4, mkChar("N"));

        result->seq = INTEGER(oseq);
    } else
        SET_VECTOR_ELT(alloc, 2, R_NilValue);

    if (param->what & WHAT_QUAL) {
        oqual = Rf_alloc3DArray(INTSXP, QUAL_LEVELS, param->n_files,
                                param->yieldSize);
        memset(INTEGER(oqual), 0, sizeof(int) * Rf_length(oqual));
        SET_VECTOR_ELT(alloc, 3, oqual);

        dimnms = NEW_LIST(3);
        Rf_setAttrib(oqual, R_DimNamesSymbol, dimnms);

        dim_elt = NEW_CHARACTER(QUAL_LEVELS);
        SET_VECTOR_ELT(dimnms, 0, dim_elt);
        SET_VECTOR_ELT(dimnms, 1, param->names);
        SET_VECTOR_ELT(dimnms, 2, R_NilValue);
        for (i = 0; i < QUAL_LEVELS; ++i) {
            qualbuf[0] = (char) (i + 33);
            SET_STRING_ELT(dim_elt, i, mkChar(qualbuf));
        }

        result->qual = INTEGER(oqual);
    } else
        SET_VECTOR_ELT(alloc, 3, R_NilValue);

    UNPROTECT(2);
    return alloc;
}

static void _mplp_setup_bam(const PILEUP_PARAM_T * param, const SPACE_T * spc,
                            PILEUP_ITER_T * plp_iter)
{
#ifdef PILEUPBAM_DEBUG
    REprintf("_mplp_setup_bam\n");
#endif
    BAM_ITER_T **mfile = plp_iter->mfile;

    for (int j = 0; j < plp_iter->n_files; ++j) {
        /* set iterator, get pileup */
        int32_t tid = bam_get_tid(mfile[j]->bfile->file->header, spc->chr);
        if (tid < 0)
            Rf_error("'%s' not in bam file %d", spc->chr, j + 1);
        mfile[j]->iter = bam_iter_query(mfile[j]->bfile->index, tid,
                                        spc->start - 1, spc->end);
    }
    plp_iter->mplp_iter =
        bam_mplp_init(plp_iter->n_files, _mplp_read_bam, (void **) mfile);
    bam_mplp_set_maxcnt(plp_iter->mplp_iter, param->max_depth);
}

static void _mplp_teardown_bam(PILEUP_ITER_T * iter)
{
#ifdef PILEUPBAM_DEBUG
    REprintf("_mplp_teardown_bam\n");
#endif
    int j;

    bam_mplp_destroy(iter->mplp_iter);
    for (j = 0; j < iter->n_files; ++j)
        bam_iter_destroy(iter->mfile[j]->iter);
}

static int _bam1(const PILEUP_PARAM_T * param, const SPACE_T * spc,
                 PILEUP_ITER_T * plp_iter, PILEUP_RESULT_T * result)
{
#ifdef PILEUPBAM_DEBUG
    REprintf("_bam1\n");
#endif
    /* A, C, G, T, N only */
    static const int nuc[] = {
        /*  A  C      G              T                          N */
        -1, 0, 1, -1, 2, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, 4
    };

    const int n_files = plp_iter->n_files, start = spc->start, end = spc->end;
    int *n_plp = plp_iter->n_plp, pos, i, j, idx = 0;

    int *opos = result->pos + result->i_yld,
        *oseq = result->seq + SEQ_LEVELS * n_files * result->i_yld,
        *oqual = result->qual + QUAL_LEVELS * n_files * result->i_yld;

    const bam_pileup1_t **plp = plp_iter->plp;
    bam_mplp_t mplp_iter = plp_iter->mplp_iter;
    int32_t tid;

    int *s0 = NULL, *q0 = NULL;

    if (param->yieldAll)
        for (i = 0; i < param->yieldSize && i < end - start + 1; ++i)
            opos[i] = start + i;

    while (param->yieldSize > idx &&
           0 < bam_mplp_auto(mplp_iter, &tid, &pos, n_plp, plp)) {
        pos += 1;
        if (pos < start || pos > end)
            continue;
        if (param->yieldAll) {
            idx = pos - start;
            if (idx >= param->yieldSize)
                break;
        } else {
            int empty = TRUE;
            for (i = 0; empty && i < n_files; ++i)
                for (j = 0; empty && j < n_plp[i]; ++j) {	/* each read */
                    const bam_pileup1_t *p = plp[i] + j;
                    if (!p->is_del || !p->is_refskip)
                        empty = FALSE;
                }
            if (empty)
                continue;
        }

        int cvg_depth = 0L;
        for (i = 0; i < n_files; ++i)
            cvg_depth += n_plp[i];
        if (param->min_depth > cvg_depth)
            continue;

        if (param->what & WHAT_SEQ)
            s0 = oseq + SEQ_LEVELS * n_files * idx;
        if (param->what & WHAT_QUAL)
            q0 = oqual + QUAL_LEVELS * n_files * idx;

        for (i = 0; i < n_files; ++i) {
            for (j = 0; j < n_plp[i]; ++j) {	/* each read */
                const bam_pileup1_t *p = plp[i] + j;
                /* filter */
                if (p->is_del || p->is_refskip)
                    continue;
                const uint8_t q = bam1_qual(p->b)[p->qpos];
                if (param->min_base_quality > q)
                    continue;
                /* query, e.g., ... */
                if (param->what & WHAT_SEQ) {
                    const int s = nuc[bam1_seqi(bam1_seq(p->b), p->qpos)];
                    if (s < 0)
                        Rf_error("unexpected nucleotide code '%d'",
                                 bam1_seqi(bam1_seq(p->b), p->qpos));
                    s0[SEQ_LEVELS * i + s] += 1;
                }
                if (param->what & WHAT_QUAL) {
                    if (QUAL_LEVELS <= q)
                        Rf_error("unexpected quality score '%ud'", q);
                    q0[QUAL_LEVELS * i + q] += 1;
                }
            }
        }
        if (!param->yieldAll)
            opos[idx] = pos;
        idx += 1;
    }
    result->i_yld += idx;
    return idx;
}

static SEXP _resize_3D_dim3(SEXP s, int n)
{
    SEXP t, dim, dimnms;
    dim = Rf_getAttrib(s, R_DimSymbol);
    dimnms = Rf_getAttrib(s, R_DimNamesSymbol);
    t = PROTECT(Rf_lengthgets(s, INTEGER(dim)[0] * INTEGER(dim)[1] * n));
    INTEGER(dim)[2] = n;
    Rf_setAttrib(t, R_DimSymbol, dim);
    Rf_setAttrib(t, R_DimNamesSymbol, dimnms);
    UNPROTECT(1);
    return t;
}

static SEXP _resize(SEXP r, int n)
{
#ifdef PILEUPBAM_DEBUG
    REprintf("_resize\n");
#endif
    SEXP s, nm = Rf_getAttrib(r, R_NamesSymbol);
    int i = 2;

    s = VECTOR_ELT(r, 1);       /* pos */
    SET_VECTOR_ELT(r, 1, Rf_lengthgets(s, n));

    s = VECTOR_ELT(r, 2);       /* seq array -- SEQ_LEVELS x n_files x n */
    if (R_NilValue != s) {
        SET_VECTOR_ELT(r, i, _resize_3D_dim3(s, n));
        SET_STRING_ELT(nm, i, STRING_ELT(nm, 2));
        ++i;
    }

    s = VECTOR_ELT(r, 3);       /* qual array -- QUAL_LEVELS x n_files x n */
    if (R_NilValue != s) {
        SET_VECTOR_ELT(r, i, _resize_3D_dim3(s, n));
        SET_STRING_ELT(nm, i, STRING_ELT(nm, 3));
        ++i;
    }

    return Rf_lengthgets(r, i);
}

static SEXP _call1(SEXP r, SEXP call)
{
#ifdef PILEUPBAM_DEBUG
    REprintf("_call1\n");
#endif
    SETCADR(call, r);           /* let's say this doesn't allocate */
    return Rf_eval(call, R_GlobalEnv);
}

static SEXP _yield1_byrange(PILEUP_PARAM_T * param, SPACE_ITER_T * spc_iter,
                            PILEUP_ITER_T * plp_iter)
{
#ifdef PILEUPBAM_DEBUG
    REprintf("_yield1_byrange\n");
#endif
    SPACE_T *spc;
    PILEUP_RESULT_T plp_result;
    SEXP res = R_NilValue, rle;
    int n_rec;

    if (NULL != (spc = _space_iter_next(spc_iter))) {
        param->yieldSize = spc->end - spc->start + 1;
        res = PROTECT(_mplp_setup_R(param, &plp_result));

        _mplp_setup_bam(param, spc, plp_iter);
        n_rec = _bam1(param, spc, plp_iter, &plp_result);
        if (param->yieldAll)
            n_rec = param->yieldSize;
        _mplp_teardown_bam(plp_iter);

        rle = _seq_rle(&n_rec, &spc->chr, 1);
        SET_VECTOR_ELT(res, 0, rle);

        res = _resize(res, n_rec);

        UNPROTECT(1);
    }

    return res;
}

static SEXP _yieldby_range(PILEUP_PARAM_T * param, SPACE_ITER_T * spc_iter,
                           PILEUP_ITER_T * plp_iter, SEXP call)
{
    SEXP result, res;
    int i;

    result = PROTECT(NEW_LIST(spc_iter->n_spc));

    for (i = 0; i < spc_iter->n_spc; ++i) {
        res = PROTECT(_yield1_byrange(param, spc_iter, plp_iter));
        if (R_NilValue == res)
            Rf_error("internal: 'spc_iter' did not yield");
        SET_VECTOR_ELT(result, i, _call1(res, call));
        UNPROTECT(1);
    }

    UNPROTECT(1);
    return result;
}

static SEXP _yield1_byposition(PILEUP_PARAM_T * param, SPACE_ITER_T * spc_iter,
                               PILEUP_ITER_T * plp_iter)
{
#ifdef PILEUPBAM_DEBUG
    REprintf("_yield1_byposition\n");
#endif
    const int yieldSize = param->yieldSize, n_spc = spc_iter->n_spc;

    SPACE_T *spc;
    PILEUP_RESULT_T plp_result;
    SEXP res = R_NilValue, rle;
    int *cnt, start, start_spc, i_spc, n_rec, i_yld = 0;

    if (NULL == (spc = _space_iter_pop(spc_iter))) {
        if (NULL == (spc = _space_iter_next(spc_iter)))
            return res;         /* early exit */
        _mplp_setup_bam(param, spc, plp_iter);
    }

    if ((spc_iter->n_spc - 1 == spc->i_spc) &&
        (param->yieldSize > spc->end - spc->start + 1))
        param->yieldSize = spc->end - spc->start + 1;

    start_spc = spc->i_spc;
    i_spc = 1;
    res = PROTECT(_mplp_setup_R(param, &plp_result));
    cnt = (int *) Calloc(n_spc, int);
    memset(cnt, 0, sizeof(int) * n_spc);

    while (spc && yieldSize > i_yld) {
        n_rec = _bam1(param, spc, plp_iter, &plp_result);
        if (param->yieldAll) {
            const int spc_width = spc->end - spc->start + 1L;
            n_rec = spc_width < param->yieldSize ? spc_width : param->yieldSize;
        }
        param->yieldSize -= n_rec;
        i_yld += n_rec;
        cnt[spc->i_spc] = i_yld;

        if (yieldSize > i_yld) {	/* next space */
            if ((spc = _space_iter_next(spc_iter))) {
                _mplp_teardown_bam(plp_iter);
                _mplp_setup_bam(param, spc, plp_iter);
                i_spc++;
            }
        }
    }

    if (i_yld) {
        rle = _seq_rle(cnt + start_spc, spc_iter->chr + start_spc, i_spc);
        SET_VECTOR_ELT(res, 0, rle);
    }
    if (spc) {
        if (i_yld) {
            start = INTEGER(VECTOR_ELT(res, 1))[i_yld - 1] + 1;
            if (spc->end >= start) {
                spc->start = start;
                _space_iter_stash(spc_iter, spc);
            }
        }
    }

    param->yieldSize = yieldSize;
    res = _resize(res, i_yld);
    Free(cnt);
    UNPROTECT(1);

    return res;
}

static SEXP _yieldby_position(PILEUP_PARAM_T * param, SPACE_ITER_T * spc_iter,
                              PILEUP_ITER_T * plp_iter, SEXP call)
{
    const int GROW_BY_ELTS = 10;
    SEXP result = R_NilValue, res;
    int i_res = 0, len, pidx;

    PROTECT_WITH_INDEX(result = NEW_LIST(0), &pidx);
    while (R_NilValue != (res = _yield1_byposition(param, spc_iter, plp_iter))) {
        PROTECT(res);
        if (Rf_length(result) == i_res) {	/* more space for results */
            len = Rf_length(result) + GROW_BY_ELTS;
            result = Rf_lengthgets(result, len);
            REPROTECT(result, pidx);
        }
        SET_VECTOR_ELT(result, i_res++, _call1(res, call));
        UNPROTECT(1);
    }
    _mplp_teardown_bam(plp_iter);	/* from _yield1_byposition */

    result = Rf_lengthgets(result, i_res);
    UNPROTECT(1);

    return result;
}

SEXP apply_pileups(SEXP files, SEXP names, SEXP space, SEXP param,
                   SEXP callback)
{
    int i;
    PILEUP_PARAM_T p;
    SPACE_ITER_T *spc_iter;
    PILEUP_ITER_T *plp_iter;
    SEXP call, result;

    if (!IS_LIST(files))
        Rf_error("'files' must be list() of BamFiles");

    p.n_files = Rf_length(files);
    p.names = names;
    for (i = 0; i < p.n_files; ++i)
        _check_isbamfile(VECTOR_ELT(files, i), "pileup");
    if (R_NilValue == space)
        Rf_error("'NULL' space not (yet) supported");
    _scan_checkparams(space, R_NilValue, R_NilValue);
    if (!Rf_isFunction(callback) || 1L != Rf_length(FORMALS(callback)))
        Rf_error("'callback' mst be a function of 1 argument");
    call = PROTECT(Rf_lang2(callback, R_NilValue));

    /* param */
    spc_iter = _space_iter_init(space);

    p.keep_flag[0] = INTEGER(_lst_elt(param, "flag", "param"))[0];
    p.keep_flag[1] = INTEGER(_lst_elt(param, "flag", "param"))[1];
    p.min_depth = INTEGER(_lst_elt(param, "minDepth", "param"))[0];
    p.max_depth = INTEGER(_lst_elt(param, "maxDepth", "param"))[0];
    p.min_base_quality = INTEGER(_lst_elt(param, "minBaseQuality", "param"))[0];
    p.min_map_quality = INTEGER(_lst_elt(param, "minMapQuality", "param"))[0];

    p.yieldSize = INTEGER(_lst_elt(param, "yieldSize", "param"))[0];
    const char *yieldBy =
        CHAR(STRING_ELT(_lst_elt(param, "yieldBy", "param"), 0));
    p.yieldBy =
        0 == strcmp(yieldBy, "range") ? YIELDBY_RANGE : YIELDBY_POSITION;
    p.yieldAll = LOGICAL(_lst_elt(param, "yieldAll", "param"))[0];

    int *what = LOGICAL(_lst_elt(param, "what", "param"));
    p.what = 0;
    if (what[0])
        p.what |= WHAT_SEQ;
    if (what[1])
        p.what |= WHAT_QUAL;

    /* data -- validate */
    plp_iter = _iter_init(files, &p);

    /* result */
    result = R_NilValue;
    if (R_NilValue == space) {  /* all */
        /* FIXME: allocate seq, but this is too big! */
        /* _bam1(n, -1, -1, max_depth,  */
        /*        _mplp_read_bam, (void **) mfile,  */
        /*        seq); */
    } else {                    /* some */
        if (YIELDBY_RANGE == p.yieldBy)
            result = _yieldby_range(&p, spc_iter, plp_iter, call);
        else
            result = _yieldby_position(&p, spc_iter, plp_iter, call);
    }

    _iter_destroy(plp_iter);
    _space_iter_destroy(spc_iter);
    UNPROTECT(1);

    return result;
}
