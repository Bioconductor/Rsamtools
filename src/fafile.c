#include <stdlib.h>             /* free */
#include "IRanges_interface.h"
#include "XVector_interface.h"
#include "fafile.h"
#include "utilities.h"

static SEXP FAFILE_TAG = NULL;

static faidx_t *_fa_tryopen(const char *fname, const char *iname)
{
    return fai_load3(fname, iname, NULL, FAI_CREATE);
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

SEXP fafile_open(SEXP filename, SEXP indexname)
{
    if (!IS_CHARACTER(filename) || 1 != Rf_length(filename))
        Rf_error("'file' must be character(1)");
    if (!IS_CHARACTER(indexname) || 1 != Rf_length(indexname))
        Rf_error("'index' must be character(1)");

    _FA_FILE *ffile = Calloc(1, _FA_FILE);
    const char
        *cfile = translateChar(STRING_ELT(filename, 0)),
        *ifile = translateChar(STRING_ELT(indexname, 0));

    ffile->index = _fa_tryopen(cfile, ifile);
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
    _checkext(ext, FAFILE_TAG, "close");
    _fafile_close(ext);
    return ext;
}

SEXP fafile_isopen(SEXP ext)
{
    SEXP ans = ScalarLogical(FALSE);
    if (NULL != FAFILE(ext)) {
        _checkext(ext, FAFILE_TAG, "isOpen");
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
    _checkext(ext, FAFILE_TAG, "countFa");
    faidx_t *fai = FAFILE(ext)->index;
    if (NULL == fai)
        Rf_error("'index' not available");
    return ScalarInteger(faidx_nseq(fai));
}

/*
 * Copy/pasted from Biostrings/src/XStringSet_io.c
 * TODO: Put it in XVector C API for sharing.
 */
static int translate(Chars_holder *seq_data, const int *lkup, int lkup_length)
{
    char *dest;
    int nbinvalid, i, j, key, val;

    /* seq_data->ptr is a const char * so we need to cast it to
       char * before we can write to it */
    dest = (char *) seq_data->ptr;
    nbinvalid = j = 0;
    for (i = 0; i < seq_data->length; i++) {
        key = (unsigned char) seq_data->ptr[i];
        if (key >= lkup_length || (val = lkup[key]) == NA_INTEGER) {
            nbinvalid++;
            continue;
        }
        dest[j++] = val;
    }
    seq_data->length = j;
    return nbinvalid;
}

SEXP scan_fa(SEXP ext, SEXP seq, SEXP start, SEXP end, SEXP type, SEXP lkup)
{
    _checkext(ext, FAFILE_TAG, "isOpen");
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

    SEXP width = PROTECT(NEW_INTEGER(n));
    int *startp = INTEGER(start), *endp = INTEGER(end),
        *widthp = INTEGER(width);
    for (int i = 0; i < n; ++i)
        widthp[i] = endp[i] - startp[i] + 1;
    const char *baseclass = CHAR(STRING_ELT(type, 0));
    char classname[13];
    snprintf(classname, sizeof(classname), "%sSet", baseclass);
    SEXP ans = PROTECT(alloc_XRawList(classname, baseclass, width));
    XVectorList_holder ans_holder = hold_XVectorList(ans);

    for (int i = 0; i < n; ++i) {
        Chars_holder ans_elt_holder = get_elt_from_XRawList_holder(&ans_holder, i);
        int len = faidx_fetch_seq2(fai, CHAR(STRING_ELT(seq, i)),
                                   startp[i] - 1, endp[i] - 1,
                                   (char *) ans_elt_holder.ptr);
        if (len == -1)
            Rf_error(" record %d (%s:%d-%d) failed", i + 1,
                     (char *) CHAR(STRING_ELT(seq, i)), startp[i], endp[i]);
        //printf("%d\n", len);
        if (len < widthp[i])
            Rf_error(" record %d (%s:%d-%d) was truncated", i + 1,
                     (char *) CHAR(STRING_ELT(seq, i)), startp[i], endp[i]);
        if (lkup != R_NilValue &&
            translate(&ans_elt_holder, INTEGER(lkup), LENGTH(lkup)) != 0)
            Rf_error(" record %d (%s:%d-%d) contains invalid DNA letters",
                     i + 1,
                     (char *) CHAR(STRING_ELT(seq, i)), startp[i], endp[i]);
    }
    UNPROTECT(2);
    return ans;
}

