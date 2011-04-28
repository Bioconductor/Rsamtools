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
    /* samtools mpileup */
    int n_files, *n_plp;
    MPLP_FILE_T **mfile;
    const bam_pileup1_t **plp;
    bam_mplp_t iter;
    /* query */
    const char *chr;
    int start, end;
    /* filter */
    int min_base_quality, min_map_quality, min_depth, max_depth;
    uint32_t keep_flag[2];
    /* yield */
    int yieldSize, yieldAll;
    YIELDBY yieldBy;
} MPLP_PARAM_T;

SEXP
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

void
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

SEXP
_mplp_setup_R(SEXP space, const int i, MPLP_PARAM_T *mparam)
{
    const int 
        start = INTEGER(VECTOR_ELT(space, 1))[i],
        end = INTEGER(VECTOR_ELT(space, 2))[i];

    mparam->chr = CHAR(STRING_ELT(VECTOR_ELT(space, 0), i));
    mparam->start = start;
    mparam->end = end;

    /* R */
    SEXP alloc = PROTECT(NEW_LIST(2)), opos, oseq;

    if (YIELDBY_RANGE == mparam->yieldBy)
        mparam->yieldSize = mparam->end - mparam->start + 1;

    opos = NEW_INTEGER(mparam->yieldSize);
    SET_VECTOR_ELT(alloc, 0, opos);
    memset(INTEGER(opos), 0, sizeof(int) * Rf_length(opos));
    oseq = Rf_alloc3DArray(INTSXP, 16, mparam->n_files,
                           mparam->yieldSize);
    memset(INTEGER(oseq), 0, sizeof(int) * Rf_length(oseq));
    SET_VECTOR_ELT(alloc, 1, oseq);

    UNPROTECT(1);

    return alloc;
}

void
_mplp_setup_bam(MPLP_PARAM_T *mparam)
{
    MPLP_FILE_T **mfile = mparam->mfile;
    int j;

    for (j = 0; j < mparam->n_files; ++j) {
        /* set iterator, get pileup */
        int32_t tid = bam_get_tid(mfile[j]->bfile->file->header, 
				  mparam->chr);
        if (tid < 0)
            Rf_error("'%s' not in bam file %d", 
		     mparam->chr, j + 1);
        mfile[j]->iter =
            bam_iter_query(mfile[j]->bfile->index, tid, 
			   mparam->start, mparam->end);
    }
    mparam->iter =
        bam_mplp_init(mparam->n_files, _mplp_read_bam, 
		      (void **) mfile);
    bam_mplp_set_maxcnt(mparam->iter, mparam->max_depth);
}

void
_mplp_teardown(MPLP_PARAM_T *mparam)
{
    int j;

    bam_mplp_destroy(mparam->iter);
    for (j = 0; j < mparam->n_files; ++j)
        bam_iter_destroy(mparam->mfile[j]->iter);
}

int
_pileup_bam1(const MPLP_PARAM_T *mparam, int *opos, int *oseq)
{
    const int n_files = mparam->n_files;
    int i, j, start = mparam->start, end = mparam->end;
    const bam_pileup1_t **plp = mparam->plp;
    int *n_plp = mparam->n_plp, pos;
    int32_t tid;

    bam_mplp_t iter = mparam->iter;

    int idx = 0;
    while (0 < bam_mplp_auto(iter, &tid, &pos, n_plp, plp) &&
           mparam->yieldSize > idx) 
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
    
    return idx;
}

SEXP
_pileup_call1(SEXP result_i, int *nrec, int i, SEXP call)
{
    int n = nrec[i % 2];
    SEXP r = VECTOR_ELT(result_i, i % 2), s, t, dim;

    if (Rf_length(VECTOR_ELT(r, 0)) != n) {
        s = VECTOR_ELT(r, 0);   /* pos */
        SET_VECTOR_ELT(r, 0, Rf_lengthgets(s, n));

        s = VECTOR_ELT(r, 1);  /* seq array -- 16 x n_files x n */
        dim = Rf_getAttrib(s, R_DimSymbol);
        t = PROTECT(Rf_lengthgets(s, INTEGER(dim)[0] * INTEGER(dim)[1] * n));
        INTEGER(dim)[2] = n;
        Rf_setAttrib(t, R_DimSymbol, dim);
        SET_VECTOR_ELT(r, 1, t);
        UNPROTECT(1);
    }

    SETCADR(call, r);
    return Rf_eval(call, R_GlobalEnv);
}

SEXP
_pileup_yieldby_range(SEXP space, MPLP_PARAM_T *mparam, SEXP call)
{
    const int nspc = Rf_length(VECTOR_ELT(space, 0));
    SEXP result, result_i, tmp;
    int *opos, *oseq, n_rec[2];

    result = PROTECT(NEW_LIST(nspc));
    result_i = PROTECT(NEW_LIST(2));

    int ispc = 0;

    tmp = _mplp_setup_R(space, ispc, mparam);
    SET_VECTOR_ELT(result_i, ispc % 2, tmp);
    opos = INTEGER(VECTOR_ELT(tmp, 0));
    oseq = INTEGER(VECTOR_ELT(tmp, 1));

    _mplp_setup_bam(mparam);
    n_rec[ispc % 2] = _pileup_bam1(mparam, opos, oseq);
    _mplp_teardown(mparam);

    /* R code needs to be evalauted by the master thread; master
     * has no implicit barrier. */
#pragma omp parallel num_threads(2) private(ispc)
    {
	for (ispc = 1; ispc < nspc; ++ispc)
	{
#pragma omp master
	    {
		tmp = _mplp_setup_R(space, ispc, mparam);
                SET_VECTOR_ELT(result_i, ispc % 2, tmp);
                opos = INTEGER(VECTOR_ELT(tmp, 0));
                oseq = INTEGER(VECTOR_ELT(tmp, 1));
	    }
#pragma omp barrier

#pragma omp master
            {
                tmp = _pileup_call1(result_i, n_rec, ispc - 1, call);
                SET_VECTOR_ELT(result, ispc - 1, tmp);
            }

#pragma omp single nowait
            {           /* next space */
		_mplp_setup_bam(mparam);
                n_rec[ispc % 2] = _pileup_bam1(mparam, opos, oseq);
                _mplp_teardown(mparam);
            }

        }
    }

    /* process final space */
    tmp = _pileup_call1(result_i, n_rec, nspc - 1, call);
    SET_VECTOR_ELT(result, nspc - 1, tmp);

    UNPROTECT(2);               /* result, result_i */

    return result;
}

SEXP
_pileup_yieldby_position(SEXP space, MPLP_PARAM_T *param, SEXP call)
{
    return R_NilValue;
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
            result = _pileup_yieldby_range(space, &mparam, call);
        else
            result = _pileup_yieldby_position(space, &mparam, call);
    }

    Free(mparam.plp);
    Free(mparam.n_plp);
    UNPROTECT(1);

    return result;
}
