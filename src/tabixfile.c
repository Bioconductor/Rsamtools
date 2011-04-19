#include <errno.h>
#include "tabixfile.h"
#include "utilities.h"

static SEXP TABIXFILE_TAG = NULL;

static void
_tabixfile_close(SEXP ext)
{
    _TABIX_FILE *tfile = TABIXFILE(ext);
    if (NULL != tfile->tabix)
	ti_close(tfile->tabix);
    tfile->tabix = NULL;
    if (NULL != tfile->iter)
	ti_iter_destroy(tfile->iter);
    tfile->iter = NULL;
}

static void
_tabixfile_finalizer(SEXP ext)
{
    if (NULL == R_ExternalPtrAddr(ext))
	return;
    _tabixfile_close(ext);
    _TABIX_FILE *tfile = TABIXFILE(ext);
    Free(tfile);
    R_SetExternalPtrAddr(ext, NULL);
}

SEXP
tabixfile_init()
{
    TABIXFILE_TAG = install("TabixFile");
    return R_NilValue;
}

SEXP 
tabixfile_open(SEXP filename, SEXP indexname)
{
    if (!IS_CHARACTER(filename) || 1L != Rf_length(filename))
	Rf_error("'filename' must be character(1)");
    if (!IS_CHARACTER(indexname) || 1L != Rf_length(indexname))
	Rf_error("'indexname' must be character(1)");

    _TABIX_FILE *tfile = Calloc(1, _TABIX_FILE);
    tfile->tabix = ti_open(translateChar(STRING_ELT(filename, 0)),
			   translateChar(STRING_ELT(indexname, 0)));
    if (NULL == tfile->tabix) {
	Free(tfile);
	Rf_error("failed to open file");
    }
    tfile->iter = NULL;

    SEXP ext = 
	PROTECT(R_MakeExternalPtr(tfile, TABIXFILE_TAG, filename));
    R_RegisterCFinalizerEx(ext, _tabixfile_finalizer, TRUE);
    UNPROTECT(1);

    return ext;
}

SEXP 
tabixfile_close(SEXP ext)
{
    _scan_checkext(ext, TABIXFILE_TAG, "close");
    _tabixfile_close(ext);
    return(ext);
}

SEXP 
tabixfile_isopen(SEXP ext)
{
    SEXP ans = ScalarLogical(FALSE);
    if (NULL != TABIXFILE(ext)) {
	_scan_checkext(ext, TABIXFILE_TAG, "isOpen");
	if (TABIXFILE(ext)->tabix)
	    ans = ScalarLogical(TRUE);
    }
    return ans;
}

SEXP
bgzip_tabix(SEXP infname, SEXP outfname, SEXP overwrite)
{
    static const int BUF_SIZE = 64 * 1024;
    void *buffer;
    int infd, oflag, outfd, cnt;
    BGZF *outp;

    if (!IS_CHARACTER(infname) || 1L != Rf_length(infname))
	Rf_error("'fromFname' must be character(1)");
    if (!IS_CHARACTER(outfname) || 1L != Rf_length(outfname))
	Rf_error("'toFname' must be character(1)");
    if (!IS_LOGICAL(overwrite) || 1L != Rf_length(overwrite))
	Rf_error("'overwrite' must be logical(1)");

    infd = open(translateChar(STRING_ELT(infname, 0)), 
		O_RDONLY);
    if (0 > infd) 
	Rf_error("opening 'fromFname': %s", strerror(errno));

    oflag = O_WRONLY | O_CREAT | O_TRUNC;
    if (!LOGICAL(overwrite)[0])
	oflag |= O_EXCL;
    outfd = open(translateChar(STRING_ELT(outfname, 0)), oflag, 0666);
    if (0 > outfd)
	Rf_error("opening 'toFname': %s", strerror(errno));

    outp = bgzf_fdopen(outfd, "w");
    if (NULL == outp)
	Rf_error("opening output 'toFname'");

    buffer = R_alloc(BUF_SIZE, sizeof(void *));
    while (0 < (cnt = read(infd, buffer, BUF_SIZE)))
	if (0 > bgzf_write(outp, buffer, cnt))
	    Rf_error("writing compressed output");
    if (0 > cnt)
	Rf_error("reading compressed output: %s", strerror(errno));

    if (0 > bgzf_close(outp))
	Rf_error("closing compressed output");
    if (-1L == close(infd))
	Rf_error("closing input after compression: %s",
		 strerror(errno));

    return outfname;
}

SEXP
index_tabix(SEXP filename, SEXP format, 
	    SEXP seq, SEXP begin, SEXP end,
	    SEXP skip, SEXP comment, SEXP zeroBased)
{
    ti_conf_t conf = ti_conf_gff;

    if (!IS_CHARACTER(filename) || 1L != Rf_length(filename))
	Rf_error("'filename' must be character(1)");
    if (1L == Rf_length(format)) {
	const char *txt = CHAR(STRING_ELT(format, 0));
	if (strcmp(txt, "gff") == 0) conf = ti_conf_gff;
	else if (strcmp(txt, "bed") == 0) conf = ti_conf_bed;
	else if (strcmp(txt, "sam") == 0) conf = ti_conf_sam;
	else if (strcmp(txt, "vcf") == 0 || 
		 strcmp(txt, "vcf4") == 0) conf = ti_conf_vcf;
	else if (strcmp(txt, "psltbl") == 0) conf = ti_conf_psltbl;
	else
	    Rf_error("format '%s' unrecognized", txt);
    } else {
	if (!IS_INTEGER(seq) || 1L != Rf_length(seq))
	    Rf_error("'seq' must be integer(1)");
	conf.sc = INTEGER(seq)[0];
	if (!IS_INTEGER(begin) || 1L != Rf_length(begin))
	    Rf_error("'begin' must be integer(1)");
	conf.bc = INTEGER(begin)[0];
	if (!IS_INTEGER(end) || 1L != Rf_length(end))
	    Rf_error("'end' must be integer(1)");
	conf.ec = INTEGER(end)[0];
    }
    
    if (IS_INTEGER(skip) && 1L == Rf_length(skip))
	conf.line_skip = INTEGER(skip)[0];
    if  (IS_CHARACTER(comment) && 1L == Rf_length(comment))
	conf.meta_char = CHAR(STRING_ELT(comment, 0))[0];
    if (IS_LOGICAL(zeroBased) && 1L == Rf_length(zeroBased))
	conf.preset |= TI_FLAG_UCSC;
    
    int res = ti_index_build(translateChar(STRING_ELT(filename, 0)),
			     &conf);

    if (-1L == res)
	Rf_error("index build failed");

    return filename;
}

SEXP
seqnames_tabix(SEXP ext)
{
    _scan_checkext(ext, TABIXFILE_TAG, "scanTabix");
    tabix_t *tabix = TABIXFILE(ext)->tabix;
    if (0 != ti_lazy_index_load(tabix))
	Rf_error("'seqnamesTabix' failed to load index");

    int n;
    const char **seqnames = ti_seqname(tabix->idx, &n);
    if (n < 0)
	Rf_error("'seqnamesTabix' found <0 (!) seqnames");
    SEXP result = PROTECT(NEW_CHARACTER(n));
    for (int i = 0; i < n; ++i)
	SET_STRING_ELT(result, i, mkChar(seqnames[i]));
    free(seqnames);
    UNPROTECT(1);
    return result;
}

SEXP 
scan_tabix(SEXP ext, SEXP space, SEXP yieldSize)
{
    const double REC_SCALE = 1.4; /* scaling factor when pre-allocated
				   * result needs to grow */
    _scan_checkparams(space, R_NilValue, R_NilValue);
    if (!IS_INTEGER(yieldSize) || 1L != Rf_length(yieldSize))
	Rf_error("'yieldSize' must be integer(1)");
    _scan_checkext(ext, TABIXFILE_TAG, "scanTabix");

    tabix_t *tabix = TABIXFILE(ext)->tabix;
    if (0 != ti_lazy_index_load(tabix))
	Rf_error("'scanTabix' failed to load index");
    
    SEXP spc = VECTOR_ELT(space, 0);
    const int
	*start = INTEGER(VECTOR_ELT(space, 1)),
	*end = INTEGER(VECTOR_ELT(space, 2)),
	nspc = Rf_length(spc);
    

    SEXP result = PROTECT(NEW_LIST(nspc));

    int buflen = 4096;
    char *buf = Calloc(buflen, char);

    for (int ispc = 0; ispc < nspc; ++ispc) {
	int totrec = INTEGER(yieldSize)[0];
	SEXP records = NEW_CHARACTER(totrec);
	SET_VECTOR_ELT(result, ispc, records); /* protect */

	int tid;
	const char *s = CHAR(STRING_ELT(spc, ispc));
	if (0 > (tid = ti_get_tid(tabix->idx, s)))
	    Rf_error("'%s' not present in tabix index", s);
	ti_iter_t iter = 
	    ti_iter_query(tabix->idx, tid, start[ispc], end[ispc]);

	int linelen;
	const char *line;
	int irec = 0;
	while (NULL != (line = ti_read(tabix, iter, &linelen))) {
	    if (totrec < irec) { /* grow */
		totrec *= REC_SCALE;
		records = Rf_lengthgets(records, totrec);
		SET_VECTOR_ELT(result, ispc, records);
	    }
	    if (linelen + 1 > buflen) {
		Free(buf);
		buflen = 2 * linelen;
		buf = Calloc(buflen, char);
	    }
	    memcpy(buf, line, linelen);
	    buf[linelen] = '\0';
	    SET_STRING_ELT(records, irec, mkChar(buf));
	    irec += 1;
	}

	ti_iter_destroy(iter);
	records = Rf_lengthgets(records, irec);
	SET_VECTOR_ELT(result, ispc, records);
    }

    Free(buf);
    UNPROTECT(1);
    return result;
}

SEXP 
yield_tabix(SEXP ext, SEXP yieldSize)
{
    if (!IS_INTEGER(yieldSize) || 1L != Rf_length(yieldSize))
	Rf_error("'yieldSize' must be integer(1)");
    _scan_checkext(ext, TABIXFILE_TAG, "scanTabix");

    tabix_t *tabix = TABIXFILE(ext)->tabix;
    ti_iter_t iter = TABIXFILE(ext)->iter;

    if (NULL == iter) {
	if (0 != ti_lazy_index_load(tabix))
	    Rf_error("'scanTabix' failed to load index");
	iter = TABIXFILE(ext)->iter = ti_iter_first();
    }

    int buflen = 4096;
    char *buf = Calloc(buflen, char);
    int linelen;
    const char *line;

    int totrec = INTEGER(yieldSize)[0];
    SEXP result = PROTECT(NEW_CHARACTER(totrec));

    int irec = 0;
    while (irec < totrec) {
	line = ti_read(tabix, iter, &linelen);
	if (NULL == line) break;
	if (linelen + 1 > buflen) {
	    Free(buf);
	    buflen = 2 * linelen;
	    buf = Calloc(buflen, char);
	}
	memcpy(buf, line, linelen);
	buf[linelen] = '\0';
	SET_STRING_ELT(result, irec, mkChar(buf));
	irec += 1;
    }

    Free(buf);
    if (irec != totrec)
	result = Rf_lengthgets(result, irec);
    UNPROTECT(1);
    return result;
}
