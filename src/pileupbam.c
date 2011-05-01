#include "samtools/bam.h"
#include "samtools/khash.h"
#include "pileupbam.h"
#include "bamfile.h"
#include "utilities.h"

typedef enum {
    YIELDBY_RANGE=0, YIELDBY_POSITION
} YIELDBY;

typedef struct {
    _BAM_FILE *bfile;
    bamFile fp;
    bam_iter_t iter;
} MPLP_FILE_T;

typedef struct {
    int i_yld;
    int *pos, *seq;
} MPLP_RESULT_T;

typedef struct {
    /* samtools mpileup */
    int n_files, *n_plp;
    MPLP_FILE_T **mfile;
    const bam_pileup1_t **plp;
    bam_mplp_t iter;
    /* query */
    int n_spc, i_spc;
    const char **chr;
    int *start, *end;
    /* filter */
    int min_base_quality, min_map_quality, min_depth, max_depth;
    uint32_t keep_flag[2];
    /* yield */
    int yieldSize, yieldAll;
    YIELDBY yieldBy;
} MPLP_PARAM_T;

static SEXP
_get_lst_elt(SEXP lst, const char *name, const char *lst_name)
{
    SEXP nms = GET_NAMES(lst);
    SEXP elt_nm = PROTECT(mkChar(name));
    int i;
    for (i = 0; i < Rf_length(nms); ++i)
        if (elt_nm == STRING_ELT(nms, i))
            break;
    UNPROTECT(1);
    if (i == Rf_length(nms))
        Rf_error("'%s' does not contain element '%s'",
                 lst_name, name);
    return VECTOR_ELT(lst, i);
}

/* from bam_aux.c; should really be exported part of library? */
KHASH_MAP_INIT_STR(s, int)

static void
_bam_init_header_hash(bam_header_t *header)
{
    if (header->hash == 0) {
        int ret, i;
        khiter_t iter;
        khash_t(s) *h;
        header->hash = h = kh_init(s);
        for (i = 0; i < header->n_targets; ++i) {
            iter = kh_put(s, h, header->target_name[i], &ret);
            kh_value(h, iter) = i;
        }
    }
}

static int
_mplp_read_bam(void *data, bam1_t *b)
{
    /* task: read one alignemnt */
    MPLP_FILE_T *mdata = (MPLP_FILE_T *) data;
    return mdata->iter ?
        bam_iter_read(mdata->fp, mdata->iter, b) :
        bam_read1(mdata->fp, b);
}

static SEXP
_mplp_setup_R(const MPLP_PARAM_T *mparam, MPLP_RESULT_T *result)
{
    SEXP alloc = PROTECT(NEW_LIST(3)), 
	nms = PROTECT(NEW_CHARACTER(3)), 
	opos, oseq;
    
    SET_STRING_ELT(nms, 0, mkChar("seqnames"));
    SET_STRING_ELT(nms, 1, mkChar("pos"));
    SET_STRING_ELT(nms, 2, mkChar("seq"));
    Rf_setAttrib(alloc, R_NamesSymbol, nms);

    /* elt 0: seqnames; handled by caller */
    opos = NEW_INTEGER(mparam->yieldSize);
    memset(INTEGER(opos), 0, sizeof(int) * Rf_length(opos));
    SET_VECTOR_ELT(alloc, 1, opos);

    oseq = Rf_alloc3DArray(INTSXP, 16, mparam->n_files,
                           mparam->yieldSize);
    memset(INTEGER(oseq), 0, sizeof(int) * Rf_length(oseq));
    SET_VECTOR_ELT(alloc, 2, oseq);

    result->i_yld = 0;
    result->pos = INTEGER(opos);
    result->seq = INTEGER(oseq);

    UNPROTECT(2);

    return alloc;
}

static void
_mplp_setup_bam(MPLP_PARAM_T *mparam)
{
    MPLP_FILE_T **mfile = mparam->mfile;
    int i_spc = mparam->i_spc, n_files = mparam->n_files;
    const char *chr = mparam->chr[i_spc];

    for (int j = 0; j < n_files; ++j) {
        /* set iterator, get pileup */
        int32_t tid = bam_get_tid(mfile[j]->bfile->file->header, chr);
        if (tid < 0)
            Rf_error("'%s' not in bam file %d", chr, j + 1);
        mfile[j]->iter =
            bam_iter_query(mfile[j]->bfile->index, tid, 
			   mparam->start[i_spc], 
			   mparam->end[i_spc]);
    }
    mparam->iter =
        bam_mplp_init(n_files, _mplp_read_bam, (void **) mfile);
    bam_mplp_set_maxcnt(mparam->iter, mparam->max_depth);
}

static void
_mplp_teardown_bam(MPLP_PARAM_T *mparam)
{
    int j;

    bam_mplp_destroy(mparam->iter);
    for (j = 0; j < mparam->n_files; ++j)
        bam_iter_destroy(mparam->mfile[j]->iter);
}

static int
_pileup_bam1(const MPLP_PARAM_T *mparam, MPLP_RESULT_T *result)
{
    const int n_files = mparam->n_files;
    int i_spc = mparam->i_spc;
    int i, j, 
	start = mparam->start[i_spc], 
	end = mparam->end[i_spc];

    int *opos = result->pos + result->i_yld,
	*oseq = result->seq + result->i_yld;
    
    const bam_pileup1_t **plp = mparam->plp;
    int *n_plp = mparam->n_plp, pos;
    int32_t tid;

    bam_mplp_t iter = mparam->iter;

    int idx = 0;
    while (mparam->yieldSize > idx &&
	   0 < bam_mplp_auto(iter, &tid, &pos, n_plp, plp))
    {
        if (pos < start || pos >= end) {
            if (mparam->yieldAll) 
                idx += 1;
            continue;
        }
        int cvg_depth = 0L;
        for (i = 0; i < n_files; ++i)
            cvg_depth += n_plp[i];
        if (mparam->min_depth > cvg_depth) {
            if (mparam->yieldAll) 
                idx += 1;
            continue;
        }
        int *s0 = oseq + 16 * n_files * idx;
        for (i = 0; i < n_files; ++i) {
            for (j = 0; j < n_plp[i]; ++j) { /* each read */
                const bam_pileup1_t *p = plp[i] + j;
                /* filter */
                if (mparam->min_map_quality > p->b->core.qual ||
                    mparam->min_base_quality > bam1_qual(p->b)[p->qpos])
                    continue;
                uint32_t test_flag = 
                    (mparam->keep_flag[0] & ~p->b->core.flag) |
                    (mparam->keep_flag[1] & p->b->core.flag);
                if (~test_flag & 2047u)
                    continue;
                /* query, e.g., ...*/
                const int s = bam1_seqi(bam1_seq(p->b), p->qpos);
                s0[16 * i + s] += 1;
            }
        }
        opos[idx] = pos;
        idx += 1;
    }

    result->i_yld += idx;
    return idx;
}

static SEXP
_pileup_call1(SEXP r, int n, SEXP call)
{
    SEXP s, t, dim;

    if (Rf_length(VECTOR_ELT(r, 1)) != n) {
        s = VECTOR_ELT(r, 1);   /* pos */
        SET_VECTOR_ELT(r, 1, Rf_lengthgets(s, n));

        s = VECTOR_ELT(r, 2);  /* seq array -- 16 x n_files x n */
        dim = Rf_getAttrib(s, R_DimSymbol);
        t = PROTECT(Rf_lengthgets(s, INTEGER(dim)[0] * INTEGER(dim)[1] * n));
        INTEGER(dim)[2] = n;
        Rf_setAttrib(t, R_DimSymbol, dim);
        SET_VECTOR_ELT(r, 2, t);
        UNPROTECT(1);
    }

    SETCADR(call, r);
    return Rf_eval(call, R_GlobalEnv);
}

static SEXP
_seq_rle(int *cnt, const char **chr, int n)
{
    int i = 0, j;
    SEXP s, t;

    for (j = 1; j < n; ++j) {
	if (0 != strcmp(chr[j], chr[j-1])) {
	    i += 1;
	    chr[i] = chr[j];
	    cnt[i] = cnt[j] - cnt[i-1];
	} 
    }
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

static SEXP
_pileup_yieldby_range(MPLP_PARAM_T *mparam, SEXP call)
{
    const int n_spc = mparam->n_spc;
    MPLP_RESULT_T mresult;
    SEXP result, result_i, res, rle;
    int i_spc = 0, n_rec[2];

    result = PROTECT(NEW_LIST(n_spc));
    result_i = PROTECT(NEW_LIST(2));

    mparam->i_spc = i_spc;
    mparam->yieldSize = mparam->end[i_spc] - mparam->start[i_spc] + 1;
    res = _mplp_setup_R(mparam, &mresult);
    SET_VECTOR_ELT(result_i, i_spc % 2, res);

    _mplp_setup_bam(mparam);
    n_rec[i_spc % 2] = _pileup_bam1(mparam, &mresult);
    _mplp_teardown_bam(mparam);

    /* R code needs to be evalauted by the master thread; master
     * has no implicit barrier. */
#pragma omp parallel num_threads(2) private(i_spc)
    {
	for (i_spc = 1; i_spc < n_spc; ++i_spc)
	{
#pragma omp master
	    {
		mparam->i_spc = i_spc;
		mparam->yieldSize = 
		    mparam->end[i_spc] - mparam->start[i_spc] + 1;
		res = _mplp_setup_R(mparam, &mresult);
                SET_VECTOR_ELT(result_i, i_spc % 2, res);
	    }
#pragma omp barrier

#pragma omp master
            {
		res = VECTOR_ELT(result_i, (i_spc - 1) % 2);
		rle = _seq_rle(&n_rec[(i_spc - 1) % 2], 
			       mparam->chr + mparam->i_spc - 1, 1);
		SET_VECTOR_ELT(res, 0, rle);
		res = _pileup_call1(res, n_rec[(i_spc - 1) % 2], call);
                SET_VECTOR_ELT(result, i_spc - 1, res);
            }

#pragma omp single nowait
            {           /* next space */
		_mplp_setup_bam(mparam);
                n_rec[i_spc % 2] = _pileup_bam1(mparam, &mresult);
                _mplp_teardown_bam(mparam);
            }
        }
    }

    /* process final space */
    res = VECTOR_ELT(result_i, (n_spc - 1) % 2);
    rle = _seq_rle(&n_rec[(n_spc - 1) % 2], 
		   mparam->chr + mparam->n_spc - 1, 1);
    SET_VECTOR_ELT(res, 0, rle);
    res = _pileup_call1(res, n_rec[(n_spc - 1) % 2], call);
    SET_VECTOR_ELT(result, n_spc - 1, res);

    UNPROTECT(2);               /* result, result_i */

    return result;
}

static SEXP
_pileup_yieldby_position(MPLP_PARAM_T *mparam, SEXP call)
{
    const int GROW_BY_ELTS = 10;
    const int n_spc = mparam->n_spc;
    const int yieldSize = mparam->yieldSize;

    SEXP result, result_i, rle;
    PROTECT_INDEX pidx;
    MPLP_RESULT_T mresult;

    int *cnt, start_spc = 0, i_yld = 0, i_result = 0;

    cnt = (int *) R_alloc(sizeof(int), n_spc);
    memset(cnt, 0, sizeof(int) * n_spc);
    PROTECT_WITH_INDEX(result = NEW_LIST(0), &pidx);

    mparam->i_spc = 0;
    _mplp_setup_bam(mparam);

    while (mparam->i_spc < n_spc) {
	if (0 == i_yld) {	/* new result */

	    if (Rf_length(result) == i_result) { /* more space for results */
		int len = Rf_length(result) + GROW_BY_ELTS;
		result = Rf_lengthgets(result, len);
		REPROTECT(result, pidx);
	    }

	    start_spc = mparam->i_spc;
	    cnt[start_spc] = 0;
	    result_i = _mplp_setup_R(mparam, &mresult);
	    SET_VECTOR_ELT(result, i_result, result_i);
	}

	i_yld += _pileup_bam1(mparam, &mresult);
	cnt[mparam->i_spc] += i_yld;

	if (i_yld == yieldSize) { /* yield */
	    rle = _seq_rle(cnt + start_spc, mparam->chr + start_spc,
			   mparam->i_spc - start_spc + 1);
	    SET_VECTOR_ELT(result_i, 0, rle);
	    result_i = _pileup_call1(result_i, i_yld, call);
	    SET_VECTOR_ELT(result, i_result++, result_i);
	    mparam->yieldSize = yieldSize;
	    i_yld = 0;

	} else {		/* next space */
	    mparam->i_spc += 1;
	    if (mparam->i_spc < mparam->n_spc) {
		mparam->yieldSize -= i_yld;
		_mplp_teardown_bam(mparam);
		_mplp_setup_bam(mparam);
	    }
	}
    }

    if (i_yld != 0) {		/* final yield */
	rle = _seq_rle(cnt + start_spc, mparam->chr + start_spc,
		       mparam->n_spc - start_spc);
	SET_VECTOR_ELT(result_i, 0, rle);
	result_i = _pileup_call1(result_i, i_yld, call);
	SET_VECTOR_ELT(result, i_result++, result_i);
	mparam->yieldSize = yieldSize;
    }
    _mplp_teardown_bam(mparam);

    result = Rf_lengthgets(result, i_result);
    UNPROTECT(1);

    return result;
}

SEXP
pileup_bam(SEXP files, SEXP space, SEXP param, 
	   SEXP callback)
{
    if (!IS_LIST(files))
        Rf_error("'files' must be list() of BamFiles");

    int i;
    MPLP_PARAM_T mparam;
    mparam.n_files = Rf_length(files);
    for (i = 0; i < mparam.n_files; ++i)
        _check_isbamfile(VECTOR_ELT(files, i), "pileup");
    _scan_checkparams(space, R_NilValue, R_NilValue);
    if (!Rf_isFunction(callback) || 1L != Rf_length(FORMALS(callback)))
        Rf_error("'callback' mst be a function of 1 argument");
    SEXP call = PROTECT(Rf_lang2(callback, R_NilValue));

    /* param */
    mparam.n_spc = Rf_length(VECTOR_ELT(space, 0));
    mparam.chr = 
	(const char **) R_alloc(sizeof(const char *), mparam.n_spc);
    for (i = 0; i < mparam.n_spc; ++i)
	mparam.chr[i] = CHAR(STRING_ELT(VECTOR_ELT(space, 0), i));
    mparam.start = INTEGER(VECTOR_ELT(space, 1));
    mparam.end = INTEGER(VECTOR_ELT(space, 2));

    mparam.keep_flag[0] =
        INTEGER(_get_lst_elt(param, "flag", "param"))[0];
    mparam.keep_flag[1] =
        INTEGER(_get_lst_elt(param, "flag", "param"))[1];
    mparam.min_depth =
        INTEGER(_get_lst_elt(param, "minDepth", "param"))[0];
    mparam.max_depth =
        INTEGER(_get_lst_elt(param, "maxDepth", "param"))[0];
    mparam.min_base_quality =
        INTEGER(_get_lst_elt(param, "minBaseQuality", "param"))[0];
    mparam.min_map_quality =
        INTEGER(_get_lst_elt(param, "minMapQuality", "param"))[0];

    mparam.yieldSize =
        INTEGER(_get_lst_elt(param, "yieldSize", "param"))[0];
    const char *yieldBy =
        CHAR(STRING_ELT(_get_lst_elt(param, "yieldBy", "param"), 0));
    mparam.yieldBy =
        0 == strcmp(yieldBy, "range") ? YIELDBY_RANGE : YIELDBY_POSITION;
    mparam.yieldAll =
        LOGICAL(_get_lst_elt(param, "yieldAll", "param"))[0];

    /* data -- validate */
    mparam.mfile =
        (MPLP_FILE_T **) R_alloc(sizeof(MPLP_FILE_T *), mparam.n_files);
    for (i = 0; i < mparam.n_files; ++i) {
        mparam.mfile[i] = (MPLP_FILE_T *) R_alloc(sizeof(MPLP_FILE_T), 1L);
        mparam.mfile[i]->bfile = BAMFILE(VECTOR_ELT(files, i));
        mparam.mfile[i]->fp = mparam.mfile[i]->bfile->file->x.bam;
        /* FIXME: the header hash should be destroyed */
        _bam_init_header_hash(mparam.mfile[i]->bfile->file->header);
    }

    mparam.plp = Calloc(mparam.n_files, const bam_pileup1_t *);
    mparam.n_plp = Calloc(mparam.n_files, int);

    /* result */
    SEXP result = R_NilValue;
    if (R_NilValue == space) {  /* all */
        /* FIXME: allocate seq, but this is too big! */
        /* _pileup_bam1(n, -1, -1, max_depth,  */
        /*            _mplp_read_bam, (void **) mfile,  */
        /*            seq); */
    } else {                    /* some */
        if (YIELDBY_RANGE == mparam.yieldBy)
            result = _pileup_yieldby_range(&mparam, call);
        else
            result = _pileup_yieldby_position(&mparam, call);
    }

    Free(mparam.plp);
    Free(mparam.n_plp);
    UNPROTECT(1);

    return result;
}
