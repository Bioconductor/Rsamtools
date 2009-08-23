#include "samtools/sam.h"
#include "io_sam.h"
#include "snap.h"
#include "utilities.h"
#include "IRanges_interface.h"

typedef enum {
	OK = 0, SEQUENCE_BUFFER_ALLOCATION_ERROR=1, 
	CIGAR_BUFFER_OVERFLOW_ERROR = 2
} _BAM_PARSE_STATUS;

typedef struct {
	int BLOCKSIZE;			  /* size to grow vectors */
	int BUF_SZ, CIGAR_BUF_SZ; /* qual / seq and cigar scratch buffer */
	char *BUF, *CIGAR_BUF;	  /* string representation of CIGAR */
	_BAM_PARSE_STATUS parse_status;
	
	bam_header_t *header;
	int nrec, idx;
	uint32_t keep_flag[2], cigar_flag;
	_SNAP_T *seq, *qual, *id;
	SEXP result;
} _BAM_DATA;


static const char *TMPL_ELT_NMS[] = {
	"qname", "flag", "rname", "strand", "pos", "width", "mapq", "cigar",
	"mrnm", "mpos", "isize", "seq", "qual"
	/* "tag", "vtype", "value" */
};

static const int N_TMPL_ELTS = sizeof(TMPL_ELT_NMS) / sizeof(const char *);

enum {
	QNAME_IDX = 0, FLAG_IDX, RNAME_IDX, STRAND_IDX, POS_IDX, WIDTH_IDX,
	MAPQ_IDX, CIGAR_IDX, MRNM_IDX, MPOS_IDX, ISIZE_IDX, SEQ_IDX,
	QUAL_IDX
};

enum { CIGAR_SIMPLE = 1 };

_BAM_DATA *
_Calloc_BAM_DATA(int blocksize, int buf_sz, int cigar_buf_sz)
{
	_BAM_DATA *bd = (_BAM_DATA *) Calloc(1, _BAM_DATA);
	bd->BLOCKSIZE = blocksize;
	bd->BUF_SZ = buf_sz;
	bd->BUF = Calloc(bd->BUF_SZ, char);
	bd->CIGAR_BUF_SZ = cigar_buf_sz;
	bd->CIGAR_BUF = Calloc(bd->CIGAR_BUF_SZ, char);
	
	bd->seq = _snap_new();
	bd->qual = _snap_new();

	return bd;
}

void
_Free_BAM_DATA(_BAM_DATA *bd)
{
	Free(bd->BUF);
	Free(bd->CIGAR_BUF);
	_snap_delete(bd->seq);
	_snap_delete(bd->qual);
	Free(bd);
}

samfile_t *
_bam_tryopen(const char * fname, const char *mode)
{
	char *fn_list = 0;
	samfile_t *sfile = samopen(fname, mode, fn_list);
	if (sfile == 0)
		Rf_error("failed to open SAM/BAM file\n  file: '%s'", fname);
	if (sfile->header == 0 || sfile->header->n_targets == 0) {
		samclose(sfile);
		Rf_error("SAM/BAM header missing or empty\n  file: '%s'", 
				 fname);
	}
	return sfile;
}

bam_index_t *
_bam_tryindexload(const char *fname)
{
	bam_index_t *index = bam_index_load(fname);
	if (index == 0)
		Rf_error("failed to load BAM index\n  file: %s", fname);
	return index;
}

static char *
_bamseq(const bam1_t *bam, char *BUF)
{
	static const char key[] = {
		'\0', 'A', 'C',  '\0',  'G', '\0', '\0', '\0', 
		'T', '\0', '\0', '\0', '\0', '\0', '\0', 'N'
	};
	
	const uint32_t len = bam->core.l_qseq;
	unsigned char *seq = bam1_seq(bam);
	for (int i = 0; i < len; ++i)
		BUF[i] = key[bam1_seqi(seq, i)];
	BUF[len] = '\0';
	if ((bam1_strand(bam) == 1))
		_reverseComplement(BUF);
	return BUF;
}

static char *
_bamqual(const bam1_t *bam, char *BUF)
{
	const uint32_t len = bam->core.l_qseq;
	unsigned char *bamq = bam1_qual(bam);
	for (int i = 0; i < len; ++i)
		BUF[i] = bamq[i] + 33;
	BUF[len] = '\0';
	if ((bam1_strand(bam) == 1))
		_reverse(BUF);
	return BUF;
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

static void
_bamfile_finalizer(SEXP externalptr)
{
	if (!R_ExternalPtrAddr(externalptr))
		return;
	samclose((samfile_t *) R_ExternalPtrAddr(externalptr));
}

SEXP
scan_bam_open(SEXP fname, SEXP mode)
{
	if (!IS_CHARACTER(fname) || LENGTH(fname) != 1)
		Rf_error("'fname' must be character(1)");
	if (!IS_CHARACTER(mode) || LENGTH(mode) != 1)
		Rf_error("'mode' must be character(1)");

	const char *file = translateChar(STRING_ELT(fname, 0));
	samfile_t *sfile = 
		_bam_tryopen(file, CHAR(STRING_ELT(mode, 0)));
	if ((sfile->type & 1) != 1) {
		samclose(sfile);
		Rf_error("'fname' is not a BAM file\n  file: %s", file);
	}
	SEXP ext = PROTECT(R_MakeExternalPtr(sfile, install("BAM_file"), fname));
	R_RegisterCFinalizerEx(ext, _bamfile_finalizer, TRUE);
	UNPROTECT(1);
	return ext;
}

void
_realloc_vectors(SEXP vec, int nrec)
{
	SEXP to, from;
	for (int i = 0; i < LENGTH(vec); ++i) {
		if (i == SEQ_IDX || i == QUAL_IDX)
			continue;
		from = VECTOR_ELT(vec, i);
		switch(TYPEOF(from)) {
		case STRSXP:
			to = NEW_STRING(nrec);
			break;
		case INTSXP:
			to = NEW_INTEGER(nrec);
			break;
		default:
			to = R_NilValue;
			break;
		}
		PROTECT(to);
		if (from != R_NilValue && LENGTH(from) > 0)
			Rf_copyVector(to, from);
		SET_VECTOR_ELT(vec, i, to);
		UNPROTECT(1);			/* 'to' protected, 'from' not */
	}
}

/* parse header */

SEXP
read_bam_header(SEXP fnames, SEXP mode, SEXP verbose)
{
	if (!IS_CHARACTER(fnames))
		error("'fnames' must be character()");
	if (!IS_CHARACTER(mode) || LENGTH(mode) != 1)
		error("'mode' must be character(1)");
	if (!IS_LOGICAL(verbose) || LENGTH(verbose) != 1)
		error("'verbose' must be logical(1)");

	SEXP ans = PROTECT(NEW_LIST(LENGTH(fnames)));
	SEXP ans_names = PROTECT(NEW_CHARACTER(LENGTH(fnames)));
	for (int i = 0; i < LENGTH(fnames); ++i) {
		samfile_t *sfile = 
			_bam_tryopen(translateChar(STRING_ELT(fnames, i)), 
						 CHAR(STRING_ELT(mode, 0)));
		bam_header_t *header = sfile->header;
		int n_elts = header->n_targets + LOGICAL(verbose)[0] ? 1 : 0;
		SEXP len = NEW_INTEGER(n_elts);
		SET_VECTOR_ELT(ans, i, len); /* protect */
		SEXP nm = PROTECT(NEW_CHARACTER(n_elts));
		for (int j = 0; j < header->n_targets; ++j) {
			SET_STRING_ELT(nm, j, mkChar(header->target_name[j]));
			INTEGER(len)[j] = header->target_len[j];
		}
			char *txt = 
				(char *) R_alloc(header->l_text + 1, sizeof(char));
			strncpy(txt, header->text, header->l_text);
			txt[header->l_text] = '\n';
			Rprintf("txt (%d): %s", header->l_text, txt);
		if (LOGICAL(verbose)[0]) {
			SET_VECTOR_ELT(ans, header->n_targets, mkChar(txt));
			SET_STRING_ELT(nm, header->n_targets, mkChar("header"));
		}
		setAttrib(len, R_NamesSymbol, nm);
		SET_STRING_ELT(ans_names, i, STRING_ELT(fnames, i));
		UNPROTECT(1);
		samclose(sfile);
	}
  
	setAttrib(ans, R_NamesSymbol, ans_names);
	UNPROTECT(2);
	return ans;
}

/* parse BAM records */

SEXP
scan_bam_template()
{
	SEXP tmpl = PROTECT(NEW_LIST(N_TMPL_ELTS));
	SET_VECTOR_ELT(tmpl, QNAME_IDX, NEW_CHARACTER(0));
	SET_VECTOR_ELT(tmpl, FLAG_IDX, NEW_INTEGER(0));
	SET_VECTOR_ELT(tmpl, RNAME_IDX, NEW_INTEGER(0));
	SET_VECTOR_ELT(tmpl, STRAND_IDX, NEW_INTEGER(0));
	SET_VECTOR_ELT(tmpl, POS_IDX, NEW_INTEGER(0));
	SET_VECTOR_ELT(tmpl, WIDTH_IDX, NEW_INTEGER(0));
	SET_VECTOR_ELT(tmpl, MAPQ_IDX, NEW_INTEGER(0));
	SET_VECTOR_ELT(tmpl, CIGAR_IDX, NEW_CHARACTER(0)); /* FIXME: cigar */
	SET_VECTOR_ELT(tmpl, MRNM_IDX, NEW_INTEGER(0));
	SET_VECTOR_ELT(tmpl, MPOS_IDX, NEW_INTEGER(0));
	SET_VECTOR_ELT(tmpl, ISIZE_IDX, NEW_INTEGER(0));
	SET_VECTOR_ELT(tmpl, SEQ_IDX, NEW_CHARACTER(0));
	SET_VECTOR_ELT(tmpl, QUAL_IDX, NEW_CHARACTER(0));

	SEXP names = PROTECT(NEW_CHARACTER(N_TMPL_ELTS));
	for (int i = 0; i < N_TMPL_ELTS; ++i)
		SET_STRING_ELT(names, i, mkChar(TMPL_ELT_NMS[i]));
	SET_ATTR(tmpl, R_NamesSymbol, names);
	UNPROTECT(2);
	return tmpl;
}

int
_scan_bam1(const bam1_t *bam, void *data)
{
	_BAM_DATA *bd = (_BAM_DATA *) data;

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
	if (bd->cigar_flag == CIGAR_SIMPLE)
	{
		if (!(n_cigar == 0 || 
			  (n_cigar == 1  && ((cigar[0] & BAM_CIGAR_MASK) == 0))))
			return 0;
	}
	if (bd->idx == bd->nrec) {
		/* reallocate */
		bd->nrec += bd->BLOCKSIZE;
		_realloc_vectors(bd->result, bd->nrec);
	}
	if (bam->core.l_qseq + 1 > bd->BUF_SZ) {
		bd->BUF_SZ = 1.4 * bam->core.l_qseq + 1;
		bd->BUF = Realloc(bd->BUF, bd->BUF_SZ, char);
		if (bd->BUF == NULL) {
			bd->parse_status |= SEQUENCE_BUFFER_ALLOCATION_ERROR;
			return -bd->idx;
		}
	}
	SEXP r = bd->result, s;
	int idx = bd->idx;
	
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
			INTEGER(s)[idx] = bam->core.tid + 1;
			break;
		case STRAND_IDX:
			INTEGER(s)[idx] = ((bam1_strand(bam) + 1) % 2) + 1;
			break;
		case POS_IDX:
			INTEGER(s)[idx] = bam->core.pos + 1;
			break;
		case WIDTH_IDX:
			INTEGER(s)[idx] = 
				bam_cigar2qlen(&bam->core, bam1_cigar(bam));
			break;
		case MAPQ_IDX:
			INTEGER(s)[idx] = bam->core.qual;
			break;
		case CIGAR_IDX:
			if (_bamcigar(cigar, n_cigar, 
						  bd->CIGAR_BUF, bd->CIGAR_BUF_SZ) < 0)
			{
				bd->parse_status |= CIGAR_BUFFER_OVERFLOW_ERROR;
				return -bd->idx;
			}
			SET_STRING_ELT(s, idx, mkChar(bd->CIGAR_BUF));
			break;
		case MRNM_IDX:
			INTEGER(s)[idx] = 
				bam->core.mtid >= 0 ? bam->core.mtid + 1 : NA_INTEGER;
			break;
		case MPOS_IDX:
			INTEGER(s)[idx] = 
				bam->core.mpos > 0 ? bam->core.mpos : NA_INTEGER;
			break;
		case ISIZE_IDX:
			INTEGER(s)[idx] = 
				bam->core.isize > 0 ? bam->core.isize : NA_INTEGER;
			break;
		case SEQ_IDX:
			_snap_append(bd->seq, _bamseq(bam, bd->BUF));
			break;
		case QUAL_IDX:
			_snap_append(bd->qual, _bamqual(bam, bd->BUF));
			break;
		default:
			break;
		}
	}
	bd->idx += 1;
	return 1;
}

int
_scan_bam_all(samfile_t *sfile, _BAM_DATA *bd)
{
	bam1_t *bam = bam_init1();
	int r = 0;
	while ((r = samread(sfile, bam)) >= 0) {
		int result = _scan_bam1(bam, bd);
		if (result < 0) return result;
	}
	return bd->idx;
}

int
_scan_bam_fetch(samfile_t *sfile, _BAM_DATA *bd, const char *fname,
				const char *space, int start, int end)
{
	int tid;
	for (tid = 0; tid < sfile->header->n_targets; ++tid) {
		if (strcmp(space, sfile->header->target_name[tid]) == 0)
			break;
	}
	if (tid == sfile->header->n_targets) {
		Rf_warning("'space' not in BAM header\n  file: %s\n  space: %s",
				   fname, space);
		return -1;
	}

	bam_index_t *bindex = _bam_tryindexload(fname);
	if (bindex == 0) {
		Rf_warning("failed to read BAM index\n  base file: %s", 
				   fname);
		return -1;
	}

	bam_fetch(sfile->x.bam, bindex, tid, start, end, bd, _scan_bam1);

	return bd->idx;
}

SEXP
_as_XStringSet(SEXP str, const char *type)
{
	SEXP nmspc = PROTECT(_get_namespace("ShortRead"));
	SEXP res, s, t;
	NEW_CALL(s, t, type, nmspc, 2);
	CSET_CDR(t, "x", str);
	CEVAL_TO(s, nmspc, res);
	UNPROTECT(1);
	return res;
}

_BAM_DATA *
_scan_bam_init(SEXP bfile, SEXP result_list, 
			   SEXP flag, SEXP isSimpleCigar)
{
	samfile_t *sfile = (samfile_t *) R_ExternalPtrAddr(bfile);
	_BAM_DATA *bdata = _Calloc_BAM_DATA(1048576, 1024, 128);
	bdata->parse_status = 0;
	bdata->header= sfile->header;
	bdata->nrec = bdata->idx = 0; 
	bdata->keep_flag[0] = INTEGER(flag)[0];
	bdata->keep_flag[1] = INTEGER(flag)[1];
	bdata->cigar_flag = LOGICAL(isSimpleCigar)[0];
	bdata->result = result_list;
	return bdata;
}

void
_scan_bam_finish(_BAM_DATA *bdata)
{
	_realloc_vectors(bdata->result, bdata->idx); /* trim */
	SEXP s;
	if ((s = VECTOR_ELT(bdata->result, STRAND_IDX)) != R_NilValue) {
		SEXP strand_lvls = PROTECT(_get_strand_levels());
		_as_factor_SEXP(s, strand_lvls);
		UNPROTECT(1);
	}
	if ((s = VECTOR_ELT(bdata->result, RNAME_IDX)) != R_NilValue) {
		_as_factor(s, (const char **) bdata->header->target_name,
				   bdata->header->n_targets);
	}
	if ((s = VECTOR_ELT(bdata->result, MRNM_IDX)) != R_NilValue) {
		_as_factor(s, (const char **) bdata->header->target_name,
				   bdata->header->n_targets);
	}
	if ((s = VECTOR_ELT(bdata->result, CIGAR_IDX)) != R_NilValue) {
		SEXP ss, tt, nmspc;
		nmspc = PROTECT(_get_namespace("Rsamtools"));
		NEW_CALL(ss, tt, "Cigar", nmspc, 2);
		CSET_CDR(tt, "cigars", s);
		CEVAL_TO(ss, nmspc, s);
		PROTECT(s);
		SET_VECTOR_ELT(bdata->result, CIGAR_IDX, s);
		UNPROTECT(2);
	}
	if ((s = VECTOR_ELT(bdata->result, SEQ_IDX)) != R_NilValue)
	{
		s = PROTECT(_snap_as_XStringSet(bdata->seq, "DNAString"));
		SET_VECTOR_ELT(bdata->result, SEQ_IDX, s);
		UNPROTECT(1);
	}
	if ((s = VECTOR_ELT(bdata->result, QUAL_IDX)) != R_NilValue)
	{
		PROTECT(s = _snap_as_PhredQuality(bdata->qual));
		SET_VECTOR_ELT(bdata->result, QUAL_IDX, s);
		UNPROTECT(1);
	}
}

SEXP
scan_bam(SEXP bfile, SEXP template_list, SEXP space, 
		 SEXP keepFlags, SEXP isSimpleCigar)
{
	if (!IS_LIST(template_list) || LENGTH(template_list) != N_TMPL_ELTS)
		Rf_error("'template' must be list(%d)", N_TMPL_ELTS);
	SEXP names = GET_ATTR(template_list, R_NamesSymbol);
	if (!IS_CHARACTER(names) || LENGTH(names) != N_TMPL_ELTS)
		Rf_error("'names(template)' must be character(%d)", 
				 N_TMPL_ELTS);
	SEXP result = PROTECT(scan_bam_template());
	for (int i = 0; i < LENGTH(names); ++i) {
		if (strcmp(TMPL_ELT_NMS[i], CHAR(STRING_ELT(names, i))) != 0)
			Rf_error("'template' names do not match scan_bam_template\n'");
		if (VECTOR_ELT(template_list, i) == R_NilValue)
			SET_VECTOR_ELT(result, i, R_NilValue);
	}
	if (!IS_INTEGER(keepFlags) || LENGTH(keepFlags) != 2)
		Rf_error("'keepFlags' must be integer(2)");
	if (!IS_LOGICAL(isSimpleCigar)  || LENGTH(isSimpleCigar) != 1)
		Rf_error("'isSimpleCigar' must be logical(1)");

	_BAM_DATA *bdata = 
		_scan_bam_init(bfile, result, keepFlags, isSimpleCigar);
	samfile_t *sfile = (samfile_t *) R_ExternalPtrAddr(bfile);

	int status;
	if (space == R_NilValue) {	/* everything */
		status = _scan_bam_all(sfile, bdata);
	} else {					/* fetch */
		const char *fname = 
			translateChar(STRING_ELT(R_ExternalPtrProtected(bfile), 0));
		const char *spc = 
			translateChar(STRING_ELT(VECTOR_ELT(space, 0), 0));
		status = _scan_bam_fetch(sfile, bdata, fname, spc,
								 INTEGER(VECTOR_ELT(space, 1))[0],
								 INTEGER(VECTOR_ELT(space, 2))[0]);
	}
	if (status < 0) {
		int idx = bdata->idx;
		_Free_BAM_DATA(bdata);
		Rf_error("failed to scan BAM\n  file: %s\n  last record: %d",
				 R_ExternalPtrProtected(bfile), idx);
	}

	_scan_bam_finish(bdata);
	_Free_BAM_DATA(bdata);
	UNPROTECT(1);
	return result;
}

/* error handling, not called directly */
void
scan_bam_cleanup()
{
	_snap_cleanup_all();
}
