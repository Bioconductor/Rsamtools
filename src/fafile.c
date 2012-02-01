#include <stdlib.h>             /* free */
#include "IRanges_interface.h"
#include "fafile.h"
#include "utilities.h"
#include "razf.h"

static SEXP FAFILE_TAG = NULL;

static faidx_t *_fa_tryopen(const char *fname)
{
    return fai_load(fname);
}

static void _fa_close(faidx_t * fai)
{
    fai_destroy(fai);
}

static void _fafile_close(SEXP ext)
{
    _FA_FILE *ffile = FAFILE(ext);
    if (NULL != ffile->index)
        _fa_close(ffile->index);
    ffile->index = NULL;
}

static void _fafile_finalizer(SEXP ext)
{
    if (NULL == R_ExternalPtrAddr(ext))
        return;
    _fafile_close(ext);
    _FA_FILE *ffile = FAFILE(ext);
    Free(ffile);
    R_SetExternalPtrAddr(ext, NULL);
}

SEXP fafile_init()
{
    FAFILE_TAG = install("FaFile");
    return R_NilValue;
}

SEXP fafile_open(SEXP filename)
{
    if (!IS_CHARACTER(filename) || 1 != Rf_length(filename))
        Rf_error("'file' must be character(1)");

    _FA_FILE *ffile = Calloc(1, _FA_FILE);
    const char *cfile = translateChar(STRING_ELT(filename, 0));

    ffile->index = _fa_tryopen(cfile);
    if (NULL == ffile->index) {
        Free(ffile);
        Rf_error("'open' index failed");
    }

    SEXP ext = PROTECT(R_MakeExternalPtr(ffile, FAFILE_TAG, filename));
    R_RegisterCFinalizerEx(ext, _fafile_finalizer, TRUE);
    UNPROTECT(1);

    return ext;
}

SEXP fafile_close(SEXP ext)
{
    _scan_checkext(ext, FAFILE_TAG, "close");
    _fafile_close(ext);
    return ext;
}

SEXP fafile_isopen(SEXP ext)
{
    SEXP ans = ScalarLogical(FALSE);
    if (NULL != FAFILE(ext)) {
        _scan_checkext(ext, FAFILE_TAG, "isOpen");
        if (NULL != FAFILE(ext)->index)
            ans = ScalarLogical(TRUE);
    }
    return ans;
}

/* fa */

SEXP index_fa(SEXP filename)
{
    if (!IS_CHARACTER(filename) || 1 != Rf_length(filename))
        Rf_error("'file' must be character(1)");

    const char *cfile = translateChar(STRING_ELT(filename, 0));
    int err = fai_build(cfile);
    if (-1 == err)
        Rf_error("'indexFa' build index failed");

    return filename;
}

SEXP n_fa(SEXP ext)
{
    _scan_checkext(ext, FAFILE_TAG, "isOpen");
    faidx_t *fai = FAFILE(ext)->index;
    if (NULL == fai)
        Rf_error("'index' not available");
    return ScalarInteger(faidx_fetch_nseq(fai));
}

SEXP scan_fa(SEXP ext, SEXP seq, SEXP start, SEXP end, SEXP lkup)
{
    _scan_checkext(ext, FAFILE_TAG, "isOpen");
    if (!IS_CHARACTER(seq))
        Rf_error("'seq' must be 'character()");
    if (!IS_INTEGER(start))
        Rf_error("'start' must be 'integer()'");
    if (!IS_INTEGER(end))
        Rf_error("'end' must be 'integer()'");
    const int n = Rf_length(seq);
    if (n != Rf_length(start) || n != Rf_length(end))
        Rf_error("'seq', 'start', and 'end' must be the same length");
    faidx_t *fai = FAFILE(ext)->index;
    if (NULL == fai)
        Rf_error("'index' not available");

    CharAEAE dna = new_CharAEAE(32767, 0);

    int *startp = INTEGER(start), *endp = INTEGER(end);
    for (int i = 0; i < n; ++i) {
        int len;
        char *seqp = faidx_fetch_seq(fai, (char *) CHAR(STRING_ELT(seq, i)),
                                     startp[i] - 1, endp[i] - 1, &len);
        if (NULL == seqp)
            Rf_error(" record %d (%s:%d-%d) failed", i + 1,
                     (char *) CHAR(STRING_ELT(seq, i)), startp[i], endp[i]);
        append_string_to_CharAEAE(&dna, seqp);
        free(seqp);
    }

    return new_XRawList_from_CharAEAE("DNAStringSet", "DNAString", &dna, lkup);
}
