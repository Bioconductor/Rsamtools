#include "samtools/sam.h"
#include "io_sam.h"
#include "encode.h"
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
	samfile_t *sfile;
	bam_header_t *header;
	int nrec, idx, irange;
	uint32_t keep_flag[2], cigar_flag;

	void *extra;
} _BAM_DATA;

typedef int (*_PARSE1_FUNC)(const bam1_t *, void *);


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

SEXP _count_bam(SEXP bfile, SEXP index, SEXP mode, 
				SEXP space, SEXP keepFlags, SEXP isSimpleCigar);

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
_Free_BAM_DATA(_BAM_DATA *bd)
{
	Free(bd->BUF);
	Free(bd->CIGAR_BUF);
	Free(bd);
}

samfile_t *
_bam_tryopen(const char *fname, const char *mode)
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

uint32_t
_bamseq(const bam1_t *bam, unsigned char *BUF)
{
	static const char key[] = {
		'\0', 'A', 'C',  '\0',  'G', '\0', '\0', '\0', 
		'T', '\0', '\0', '\0', '\0', '\0', '\0', 'N'
	};
	
	uint32_t len = bam->core.l_qseq;
	unsigned char *seq = bam1_seq(bam);
	for (int i = 0; i < len; ++i)
		BUF[i] = key[bam1_seqi(seq, i)];
	if ((bam1_strand(bam) == 1))
		_reverseComplement(BUF, len);
	return len;
}

uint32_t
_bamqual(const bam1_t *bam, unsigned char *BUF)
{
	uint32_t len = bam->core.l_qseq;
	unsigned char *bamq = bam1_qual(bam);
	for (int i = 0; i < len; ++i)
		BUF[i] = bamq[i] + 33;
	if ((bam1_strand(bam) == 1))
		_reverse(BUF, len);
	return len;
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
	if (NULL == R_ExternalPtrAddr(externalptr))
		return;
	samclose((samfile_t *) R_ExternalPtrAddr(externalptr));
	R_SetExternalPtrAddr(externalptr, NULL);
}

SEXP
_scan_bam_open(SEXP fname, SEXP mode)
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
scan_bam_close(SEXP bfile)
{
	_bamfile_finalizer(bfile);
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
		int n_elts = 
			header->n_targets + (LOGICAL(verbose)[0] ? 1 : 0);
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
		txt[header->l_text] = '\0';
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

/* general functions */

void
_scan_check_params(SEXP space, SEXP keepFlags, SEXP isSimpleCigar)
{
	if (!IS_INTEGER(keepFlags) || LENGTH(keepFlags) != 2)
		Rf_error("'keepFlags' must be integer(2)");
	if (!IS_LOGICAL(isSimpleCigar)  || LENGTH(isSimpleCigar) != 1)
		Rf_error("'isSimpleCigar' must be logical(1)");
}

_BAM_DATA *
_init_BAM_DATA(SEXP bfile, SEXP index, SEXP flag, SEXP isSimpleCigar)
{
	_BAM_DATA *bdata = _Calloc_BAM_DATA(1048576, 1024, 128);
	bdata->parse_status = 0;
	bdata->sfile = (samfile_t *) R_ExternalPtrAddr(bfile);
	bdata->header= bdata->sfile->header;
	bdata->nrec = bdata->idx = bdata->irange = 0; 
	bdata->keep_flag[0] = INTEGER(flag)[0];
	bdata->keep_flag[1] = INTEGER(flag)[1];
	bdata->cigar_flag = LOGICAL(isSimpleCigar)[0];
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
	while ((r = samread(bd->sfile, bam)) >= 0) {
		int result = (*parse1)(bam, bd);
		if (result < 0) return result;
	}
	return bd->idx;
}

int
_scan_bam_fetch(_BAM_DATA *bd, const char *fname, const char *indexfname,
				SEXP space, int* start, int* end,
				_PARSE1_FUNC parse1)
{
	int tid;
	samfile_t *sfile = bd->sfile;
	int n_tot = 0;
	
	for (int irange = 0; irange < LENGTH(space); ++irange) {
		const char* spc = translateChar(STRING_ELT(space, irange));
		for (tid = 0; tid < sfile->header->n_targets; ++tid) {
			if (strcmp(spc, sfile->header->target_name[tid]) == 0)
				break;
		}
		if (tid == sfile->header->n_targets) {
			Rf_warning("'space' not in BAM header\n  file: %s\n  space: %s",
					   fname, spc);
			return -1;
		}

		bam_index_t *bindex = _bam_tryindexload(indexfname);
		if (bindex == 0) {
			Rf_warning("failed to read BAM index\n  base file: %s", 
					   fname);
			return -1;
		}

		bam_fetch(sfile->x.bam, bindex, tid, 
				  start[irange], end[irange], bd, parse1);
		n_tot += bd->idx;
		bd->irange += 1;
		bd->idx = 0;
	}
	return n_tot;
}

int
_do_scan_bam(_BAM_DATA *bdata, SEXP bfile, SEXP index,
			 SEXP space, _PARSE1_FUNC parse1)
{
	int status;

	if (space == R_NilValue) {	/* everything */
		status = _scan_bam_all(bdata, parse1);
	} else {					/* fetch */
		const char *fname = 
			translateChar(STRING_ELT(R_ExternalPtrProtected(bfile), 0));
		const char *indexfname =
			translateChar(STRING_ELT(index, 0));
		status = _scan_bam_fetch(bdata, fname, indexfname,
								 VECTOR_ELT(space, 0),
								 INTEGER(VECTOR_ELT(space, 1)),
								 INTEGER(VECTOR_ELT(space, 2)),
								 parse1);
	}
	return status;
}

/* parse */

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
_scan_bam_parse1(const bam1_t *bam, void *data)
{
	_BAM_DATA *bd = (_BAM_DATA *) data;
	if (FALSE == _bam_filter(bam, bd))
		return 0;

	SEXP s, r = VECTOR_ELT((SEXP) bd->extra, bd->irange);
	int idx = bd->idx, start, len;
	
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
			if (_bamcigar(bam1_cigar(bam), bam->core.n_cigar,
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
			start = TRUELENGTH(s);
			len = _bamseq(bam, RAW(VECTOR_ELT(s, 0)) + start);
			SET_TRUELENGTH(s, start + len);
			INTEGER(VECTOR_ELT(s, 1))[idx] = len;
			break;
		case QUAL_IDX:
			start = TRUELENGTH(s);
			len = _bamqual(bam, RAW(VECTOR_ELT(s, 0)) + start);
			SET_TRUELENGTH(s, start + len);
			INTEGER(VECTOR_ELT(s, 1))[idx] = len;
			break;
		default:
			break;
		}
	}
	bd->idx += 1;
	return 1;
}

SEXP
_scan_bam_result_init(SEXP count, SEXP template_list, SEXP names)
{
 	int nrange = LENGTH(VECTOR_ELT(count, 0));
	SEXP result = PROTECT(NEW_LIST(nrange));
	/* result: 
	   range1: tmpl1, tmpl2...
	   range2: tmpl1, tmpl2...
	   ...
	*/
	for (int irange = 0; irange < nrange; ++irange) 
	{
		int n_read = INTEGER(VECTOR_ELT(count, 0))[irange],
			n_nucs = INTEGER(VECTOR_ELT(count, 1))[irange];

		SEXP tmpl = PROTECT(scan_bam_template());
		for (int i = 0; i < LENGTH(names); ++i) {
			if (VECTOR_ELT(template_list, i) == R_NilValue)
				SET_VECTOR_ELT(tmpl, i, R_NilValue);
			else {
				if (SEQ_IDX == i || QUAL_IDX == i) {
					SEXP lst = PROTECT(NEW_LIST(2));
					SET_TRUELENGTH(lst, 0); /* nucleotide count */
					SET_VECTOR_ELT(lst, 0, NEW_RAW(n_nucs));
					SET_VECTOR_ELT(lst, 1, NEW_INTEGER(n_read));
					SET_VECTOR_ELT(tmpl, i, lst);
					UNPROTECT(1);
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
_as_XStringSet(SEXP lst, const char *baseclass)
{
	SEXP seq = VECTOR_ELT(lst, 0);
	SEXP width = VECTOR_ELT(lst, 1);

	ENCODE_FUNC encode = _encoder(baseclass);
	char *str = (char *) RAW(seq);
	for (int i = 0; i < LENGTH(seq); ++i)		
		str[i] = encode(str[i]);
	SEXP ptr = PROTECT(new_SequencePtr("RawPtr", seq));
	SEXP xstring = PROTECT(new_XSequence(baseclass, ptr, 0, LENGTH(seq)));

	SEXP start = PROTECT(NEW_INTEGER(LENGTH(width)));
	int s = 1;
	for (int i = 0; i < LENGTH(width); ++i) {
		INTEGER(start)[i] = s;
		s += INTEGER(width)[i];
	}
	SEXP irange = 
		PROTECT(new_IRanges("IRanges", start, width, R_NilValue));

	const int XSETCLASS_BUF = 40;
	if (strlen(baseclass) > XSETCLASS_BUF - 4)
		error("Rsamtools internal error: *Set buffer too small");
	char xsetclass[XSETCLASS_BUF];
	snprintf(xsetclass, XSETCLASS_BUF, "%sSet", baseclass);
	SEXP xclassdef = PROTECT(MAKE_CLASS(xsetclass));
	SEXP xstringset = PROTECT(NEW_OBJECT(xclassdef));
        
	SET_SLOT(xstringset, install("super"), xstring);
	SET_SLOT(xstringset, install("ranges"), irange);

	UNPROTECT(6);
	return xstringset;
}

SEXP
_as_PhredQuality(SEXP lst)
{
	SEXP xstringset = PROTECT(_as_XStringSet(lst, "BString"));

	SEXP s, t, nmspc, result;
	nmspc = PROTECT(_get_namespace("Rsamtools"));
	NEW_CALL(s, t, "PhredQuality", nmspc, 2);
	CSET_CDR(t, "x", xstringset);
	CEVAL_TO(s, nmspc, result);
	UNPROTECT(2);
	return result;
}

void
_scan_bam_finish1range(_BAM_DATA *bdata, SEXP result)
{
	SEXP s;
	if ((s = VECTOR_ELT(result, STRAND_IDX)) != R_NilValue) {
		SEXP strand_lvls = PROTECT(_get_strand_levels());
		_as_factor_SEXP(s, strand_lvls);
		UNPROTECT(1);
	}
	if ((s = VECTOR_ELT(result, RNAME_IDX)) != R_NilValue) {
		_as_factor(s, (const char **) bdata->header->target_name,
				   bdata->header->n_targets);
	}
	if ((s = VECTOR_ELT(result, MRNM_IDX)) != R_NilValue) {
		_as_factor(s, (const char **) bdata->header->target_name,
				   bdata->header->n_targets);
	}
	if ((s = VECTOR_ELT(result, CIGAR_IDX)) != R_NilValue) {
		SEXP ss, tt, nmspc;
		nmspc = PROTECT(_get_namespace("Rsamtools"));
		NEW_CALL(ss, tt, "Cigar", nmspc, 2);
		CSET_CDR(tt, "cigars", s);
		CEVAL_TO(ss, nmspc, s);
		PROTECT(s);
		SET_VECTOR_ELT(result, CIGAR_IDX, s);
		UNPROTECT(2);
	}
	if ((s = VECTOR_ELT(result, SEQ_IDX)) != R_NilValue)
	{
		s = PROTECT(_as_XStringSet(s, "DNAString"));
		SET_VECTOR_ELT(result, SEQ_IDX, s);
		UNPROTECT(1);
	}
	if ((s = VECTOR_ELT(result, QUAL_IDX)) != R_NilValue)
	{
		PROTECT(s = _as_PhredQuality(s));
		SET_VECTOR_ELT(result, QUAL_IDX, s);
		UNPROTECT(1);
	}
}

void
_scan_bam_finish(_BAM_DATA *bdata)
{
	SEXP result = (SEXP) bdata->extra;
	for (int irange = 0; irange < LENGTH(result); ++irange)
		_scan_bam_finish1range(bdata, VECTOR_ELT(result, irange));
}

SEXP
scan_bam(SEXP fname, SEXP index, SEXP mode, SEXP template_list, 
		 SEXP space, SEXP keepFlags, SEXP isSimpleCigar)
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
	_scan_check_params(space, keepFlags, isSimpleCigar);

	SEXP bfile = PROTECT(_scan_bam_open(fname, mode));
	SEXP count = PROTECT(_count_bam(bfile, index, mode, 
									space, keepFlags, isSimpleCigar));
	if (R_NilValue == count) {
		scan_bam_close(bfile);
		UNPROTECT(2);
		Rf_error("Failed to scan BAM\n  file: %s", 
				 translateChar(STRING_ELT(fname, 0)));
	}
	if (R_NilValue == space) {
		/* bam file invalid if fully scanned */
		scan_bam_close(bfile);
		UNPROTECT(2); PROTECT(count);
		bfile = PROTECT(_scan_bam_open(fname, mode));
	}

	SEXP result = 
		PROTECT(_scan_bam_result_init(count, template_list, names));
	_BAM_DATA *bdata = _init_BAM_DATA(bfile, index, keepFlags, isSimpleCigar);
	bdata->extra = (void *) result;

	int status = 
		_do_scan_bam(bdata, bfile, index, space, _scan_bam_parse1);
	if (status < 0) {
		int idx = bdata->idx;
		const char *fname0 = 
			translateChar(STRING_ELT(R_ExternalPtrProtected(bfile), 0));
		scan_bam_close(bfile);
		_Free_BAM_DATA(bdata);
		Rf_error("failed to scan BAM\n  file: %s\n  last record: %d",
				 fname0, idx);
	}

	_scan_bam_finish(bdata);
	_Free_BAM_DATA(bdata);
	scan_bam_close(bfile);
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
	INTEGER(VECTOR_ELT(cnt, 1))[bd->irange] += bam->core.l_qseq;
	return 1;
}

SEXP
_count_bam(SEXP bfile, SEXP index, SEXP mode, 
		   SEXP space, SEXP keepFlags, SEXP isSimpleCigar)
{
	SEXP result = PROTECT(NEW_LIST(2));
	int spc_length = 
		(R_NilValue == space) ? 1 : LENGTH(VECTOR_ELT(space, 0));
	SET_VECTOR_ELT(result, 0, NEW_INTEGER(spc_length));
	SET_VECTOR_ELT(result, 1, NEW_INTEGER(spc_length));
	for (int i = 0; i < spc_length; ++i) {
		INTEGER(VECTOR_ELT(result, 0))[i] =
			INTEGER(VECTOR_ELT(result, 1))[i] = 0;
	}

	SEXP nms = PROTECT(NEW_CHARACTER(2));
	SET_STRING_ELT(nms, 0, mkChar("records"));
	SET_STRING_ELT(nms, 1, mkChar("nucleotides"));
	setAttrib(result, R_NamesSymbol, nms);
	UNPROTECT(1);

	_BAM_DATA *bdata = _init_BAM_DATA(bfile, index, keepFlags, isSimpleCigar);
	bdata->extra = result;

	int status = 
		_do_scan_bam(bdata, bfile, index, space, _count_bam1);
	if (status < 0) 
		result = R_NilValue;

	_Free_BAM_DATA(bdata);
	UNPROTECT(1);
	return result;
}

SEXP
count_bam(SEXP fname, SEXP index, SEXP mode, 
		  SEXP space, SEXP keepFlags, SEXP isSimpleCigar)
{	
	_scan_check_params(space, keepFlags, isSimpleCigar);
	SEXP bfile = PROTECT(_scan_bam_open(fname, mode));
	SEXP count = PROTECT(_count_bam(bfile, index, mode, space, 
									keepFlags, isSimpleCigar));
	scan_bam_close(bfile);
	UNPROTECT(2);

	if (R_NilValue == count)
		Rf_error("failed to count BAM\n  file: %s", fname);
	return count;
}

void
scan_bam_cleanup()
{
	/* placeholder */
}
