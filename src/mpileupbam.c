#include "samtools/bam.h"
#include "samtools/khash.h"
#include "mpileupbam.h"
#include "bamfile.h"
#include "utilities.h"

typedef struct {
    bamFile fp;
    bam_iter_t iter;
}  MPLP_DATA_T;

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
    MPLP_DATA_T *mdata = (MPLP_DATA_T *) data;
    return mdata->iter ?
	bam_iter_read(mdata->fp, mdata->iter, b) :
	bam_read1(mdata->fp, b);
}

void
_mpileup_bam1(const int n, const int start, const int end, 
	      const int max_depth, bam_plp_auto_f fun, void **data,
	      int *seq)
{
    int *n_plp = Calloc(n, int), i, j, tid, pos;
    const bam_pileup1_t **plp = Calloc(n, const bam_pileup1_t *);
    bam_mplp_t iter = bam_mplp_init(n, fun, data);
    bam_mplp_set_maxcnt(iter, max_depth);
    while (bam_mplp_auto(iter, &tid, &pos, n_plp, plp) > 0) {
	if (pos < start || pos >= end) continue;
	for (i = 0; i < n; ++i) { /* each file */
	    int *s0 = seq +
		(i * (end - start + 1) + (pos - start)) * 16;
	    for (j = 0; j < n_plp[i]; ++j) { /* each read */
		/* query individual pileup, e.g., ...*/
		const bam_pileup1_t *p = plp[i] + j;
		const int s = bam1_seqi(bam1_seq(p->b), p->qpos);
		s0[s] += 1;
	    }
	}
    }
    bam_mplp_destroy(iter);
    Free(plp); Free(n_plp);
}

SEXP
mpileup_bam(SEXP files, SEXP space, SEXP param, 
	    SEXP callback)
{
    int i, j;

    if (!IS_LIST(files))
	Rf_error("'files' must be list() of BamFiles");
    const int n = Rf_length(files);
    for (i = 0; i < n; ++i)
	_check_isbamfile(VECTOR_ELT(files, i), "mpileup");
    _scan_checkparams(space, R_NilValue, R_NilValue);
    if (!Rf_isFunction(callback) || 1L != Rf_length(FORMALS(callback)))
	Rf_error("'callback' mst be a function of 1 argument");
    SEXP call = PROTECT(Rf_lang2(callback, R_NilValue));
    /* FIXME: check param */

    /* param */
    const int max_depth =
	INTEGER(_get_lst_elt(param, "maxDepth", "param"))[0];

    /* data -- validate */
    MPLP_DATA_T **data =
	(MPLP_DATA_T **) R_alloc(sizeof(MPLP_DATA_T *), n);
    for (i = 0; i < n; ++i) {
	_BAM_FILE *bfile = BAMFILE(VECTOR_ELT(files, i));
	data[i] = (MPLP_DATA_T *) R_alloc(sizeof(MPLP_DATA_T), 1L);
	data[i]->fp = bfile->file->x.bam;
	/* FIXME: the header hash should be destroyed */
	_bam_init_header_hash(bfile->file->header);
    }

    /* result */
    SEXP result;
    /* spaces */
    if (R_NilValue == space) {
	/* FIXME: allocate seq, but this is too big! */
	/* _mpileup_bam1(n, -1, -1, max_depth,  */
	/* 	      _mplp_read_bam, (void **) data,  */
	/* 	      seq); */
    } else {
	SEXP spc = VECTOR_ELT(space, 0);
	const int
	    *start = INTEGER(VECTOR_ELT(space, 1)),
	    *end = INTEGER(VECTOR_ELT(space, 2)),
	    nspc = Rf_length(spc);
	result = PROTECT(NEW_LIST(nspc));
	int **resi = Calloc(n, int *);
	for (i = 0; i < nspc; ++i) {
	    SEXP res1 = NEW_LIST(2); /* eventually: seq, qual */
	    SET_VECTOR_ELT(result, i, res1);

	    /* allocate results space */
	    SEXP tmp = 
		Rf_alloc3DArray(INTSXP, 16, end[i] - start[i] + 1, n);
	    memset(INTEGER(tmp), 0, sizeof(int) * Rf_length(tmp));
	    SET_VECTOR_ELT(res1, 0, tmp);

	    const char *s = CHAR(STRING_ELT(spc, i));
	    for (j = 0; j < n; ++j) {
		/* set iterator */
		_BAM_FILE *bfile  = BAMFILE(VECTOR_ELT(files, j));
		int32_t tid = bam_get_tid(bfile->file->header, s);
		if (tid < 0)
		    Rf_error("'%s' not in bam file %d", s, j + 1);
		data[j]->iter =
		    bam_iter_query(bfile->index, tid, start[i], end[i]);
	    }
	    _mpileup_bam1(n, start[i], end[i], max_depth, 
			  _mplp_read_bam, (void **) data, 
			  INTEGER(VECTOR_ELT(res1, 0)));
	    for (j = 0; j < n; ++j)
		bam_iter_destroy(data[j]->iter);

	    SETCADR(call, res1);
	    SET_VECTOR_ELT(result, i, Rf_eval(call, R_GlobalEnv));
	}

	Free(resi);
    }

    UNPROTECT(2);
    return result;
}
