#include "samtools/bam.h"
#include "samtools/khash.h"
#include "mpileupbam.h"
#include "bamfile.h"
#include "utilities.h"

typedef struct {
    _BAM_FILE *bfile;
    bamFile fp;
    bam_iter_t iter;
} MPLP_FILE_T;

typedef struct {
    MPLP_FILE_T **mfile;
    int n_files;
    const bam_pileup1_t **plp;
    int *n_plp;
    /* query */
    const char *chr;
    int start, end;
    /* filter */
    int min_base_quality, min_map_quality, min_depth, max_depth;
    uint32_t keep_flag[2];
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
_mplp_setup(SEXP space, int i, MPLP_PARAM_T *mparam)
{
    const int 
	start = INTEGER(VECTOR_ELT(space, 1))[i],
	end = INTEGER(VECTOR_ELT(space, 2))[i];
    mparam->chr = CHAR(STRING_ELT(VECTOR_ELT(space, 0), i));
    mparam->start = start;
    mparam->end = end;

    mparam->plp = Calloc(mparam->n_files, const bam_pileup1_t *);
    mparam->n_plp = Calloc(mparam->n_files, int);

    SEXP tmp = 
	Rf_alloc3DArray(INTSXP, 16, end - start + 1, mparam->n_files);
    memset(INTEGER(tmp), 0, sizeof(int) * Rf_length(tmp));
    return tmp;
}

void
_mplp_teardown(MPLP_PARAM_T *mparam)
{
    Free(mparam->plp);
    Free(mparam->n_plp);
}

void
_mpileup_bam1(const MPLP_PARAM_T *mparam, int *seq)
{
    MPLP_FILE_T **mfile = mparam->mfile;
    const int n_files = mparam->n_files;
    int i, j, start = mparam->start, end = mparam->end;
    const char *chr = mparam->chr;
    int32_t tid;

    for (j = 0; j < n_files; ++j) {
	/* set iterator, get pileup */
	tid = bam_get_tid(mfile[j]->bfile->file->header, chr);
	if (tid < 0)
	    Rf_error("'%s' not in bam file %d", chr, j + 1);
	mfile[j]->iter =
	    bam_iter_query(mfile[j]->bfile->index, tid, start, end);
    }

    const bam_pileup1_t **plp = mparam->plp;
    int *n_plp = mparam->n_plp, pos;

    bam_mplp_t iter = 
	bam_mplp_init(n_files, _mplp_read_bam, (void **) mfile);
    bam_mplp_set_maxcnt(iter, mparam->max_depth);
    while (bam_mplp_auto(iter, &tid, &pos, n_plp, plp) > 0) {
	if (pos < start || pos >= end) 
	    continue;
	int cvg_depth = 0L;
	for (i = 0; i < n_files; ++i)
	    cvg_depth += n_plp[i];
	if (mparam->min_depth > cvg_depth)
	    continue;
	for (i = 0; i < n_files; ++i) {
	    int *s0 =
		seq + (i * (end - start + 1) + (pos - start)) * 16;
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
		s0[s] += 1;
	    }
	}
    }
    bam_mplp_destroy(iter);

    for (j = 0; j < mparam->n_files; ++j)
	bam_iter_destroy(mfile[j]->iter);
}

SEXP
mpileup_bam(SEXP files, SEXP space, SEXP param, 
	    SEXP callback)
{
    if (!IS_LIST(files))
	Rf_error("'files' must be list() of BamFiles");

    int i;
    MPLP_PARAM_T mparam;
    mparam.n_files = Rf_length(files);
    for (i = 0; i < mparam.n_files; ++i)
	_check_isbamfile(VECTOR_ELT(files, i), "mpileup");
    _scan_checkparams(space, R_NilValue, R_NilValue);
    if (!Rf_isFunction(callback) || 1L != Rf_length(FORMALS(callback)))
	Rf_error("'callback' mst be a function of 1 argument");
    SEXP call = PROTECT(Rf_lang2(callback, R_NilValue));

    /* param */
    /* FIXME: check param */
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

    /* result */
    SEXP result;
    if (R_NilValue == space) {	/* all */
	/* FIXME: allocate seq, but this is too big! */
	/* _mpileup_bam1(n, -1, -1, max_depth,  */
	/* 	      _mplp_read_bam, (void **) mfile,  */
	/* 	      seq); */
    } else {			/* some */
        const int nspc = Rf_length(VECTOR_ELT(space, 0));
	result = PROTECT(NEW_LIST(nspc));

	SEXP resi[2], tmp;	/* space for results */
	int *itmp;
	for (i = 0; i < 2; ++i)
	    resi[i] = PROTECT(NEW_LIST(2));

	i = 0;

	/* R code needs to be evalauted by the master thread; master
	 * has no implicit barrier. */
#pragma omp parallel num_threads(2)
	{
#pragma omp master
	    {
		tmp = _mplp_setup(space, i, &mparam);
		SET_VECTOR_ELT(resi[ i % 2 ], 0, tmp);
		itmp = INTEGER(tmp);
		/* first space. No parallelism here; eval on master */
		_mpileup_bam1(&mparam, itmp);
	    }
#pragma omp barrier

	for (i = 1; i < nspc; ++i) {
#pragma omp master
	    {
		tmp = _mplp_setup(space, i, &mparam);
		SET_VECTOR_ELT(resi[ i % 2 ], 0, tmp);
		itmp = INTEGER(tmp);
	    }
#pragma omp barrier

#pragma omp single nowait
		{		/* next space */
		    _mpileup_bam1(&mparam, itmp);
		}
#pragma omp master
		{
		    SETCADR(call, resi[ (i-1) % 2]);
		    SET_VECTOR_ELT(result, i-1, Rf_eval(call, R_GlobalEnv));
		}
	    }
	}

	/* process final space */
	SETCADR(call, resi[ (nspc-1) % 2 ]);
	SET_VECTOR_ELT(result, nspc-1, Rf_eval(call, R_GlobalEnv));

	_mplp_teardown(&mparam);
	UNPROTECT(3);		/* resi */
    }

    UNPROTECT(1);
    return result;
}
