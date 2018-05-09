#include "Rdefines.h"
#include "bamfile.h"
#include "bambuffer.h"
#include "bam_data.h"
#include "scan_bam_data.h"
#include "utilities.h"
#include "io_sam.h"     /* _bam_check_template_list */

static SEXP BAMBUFFER_TAG = NULL;

BAM_BUFFER bambuffer_new(int n, int as_mates)
{
    BAM_BUFFER buf = Calloc(1, _BAM_BUFFER);
    buf->i = 0;
    buf->n = n;
    buf->buffer = Calloc(n, bam1_t *);
    if (as_mates) {
        buf->as_mates = TRUE;
        buf->mates = Calloc(n, int);
        buf->partition = Calloc(n, int);
    }
    return buf;
}

void bambuffer_push(BAM_BUFFER buf, const bam1_t *bam)
{
    if (buf->i == buf->n) {
        buf->n *= 1.3;
        buf->buffer = Realloc(buf->buffer, buf->n, bam1_t *);
        if (buf->as_mates) {
            buf->mates = Realloc(buf->mates, buf->n, int);
            buf->partition = Realloc(buf->partition, buf->n, int);
        }
    }
    buf->buffer[buf->i] = bam_dup1(bam);
    if (buf->as_mates) {
        buf->mates[buf->i] = buf->mate_flag;
        buf->partition[buf->i] = buf->partition_id;
    }
    buf->i += 1;
}

void _bambuffer_reset(BAM_BUFFER buf)
{
    for (int i = 0; i < buf->i; ++i)
        bam_destroy1(buf->buffer[i]);
    buf->i = 0;
}

void bambuffer_free(BAM_BUFFER buf)
{
    _bambuffer_reset(buf);
    Free(buf->buffer);
    if (buf->as_mates) {
        Free(buf->mates);
        Free(buf->partition);
    }
    Free(buf);
}

static void _bambuffer_finalizer(SEXP ext)
{
    if (NULL == R_ExternalPtrAddr(ext))
        return;
    BAM_BUFFER buf = BAMBUFFER(ext);
    bambuffer_free(buf);
    R_SetExternalPtrAddr(ext, NULL);
}

/* SEXP interface */

SEXP bambuffer_init()
{
    BAMBUFFER_TAG = install("BamBuffer");
    return R_NilValue;
}

SEXP bambuffer(int yieldSize, int as_mates)
{
    BAM_BUFFER buf = bambuffer_new(yieldSize, as_mates);
    SEXP ext = PROTECT(R_MakeExternalPtr(buf, BAMBUFFER_TAG, NULL));
    R_RegisterCFinalizerEx(ext, _bambuffer_finalizer, TRUE);
    UNPROTECT(1);
    return ext;
}

SEXP bambuffer_length(SEXP bufext)
{
    _checkext(bufext, BAMBUFFER_TAG, "bamBuffer 'parse'");
    return ScalarInteger(BAMBUFFER(bufext)->i);
}

SEXP bambuffer_parse(SEXP ext, SEXP regions, SEXP keepFlags, SEXP isSimpleCigar,
                     SEXP tagFilter, SEXP mapqFilter, SEXP bufext,
                     SEXP reverseComplement, SEXP templateList)
{
    _check_isbamfile(ext, "bamBuffer, 'parse'");
    _checkparams(regions, keepFlags, isSimpleCigar);
    _checkext(bufext, BAMBUFFER_TAG, "bamBuffer 'parse'");
    if (!(IS_LOGICAL(reverseComplement) && (1L == LENGTH(reverseComplement))))
        Rf_error("'reverseComplement' must be logical(1)");
    _bam_check_template_list(templateList);
    
    SEXP names = GET_ATTR(templateList, R_NamesSymbol);
    SEXP result =
        PROTECT(_scan_bam_result_init(templateList, names, R_NilValue,
                                      BAMFILE(ext)));
    SCAN_BAM_DATA sbd = _init_SCAN_BAM_DATA(result);
    BAM_DATA bd = _init_BAM_DATA(ext, R_NilValue, keepFlags, isSimpleCigar,
                                 tagFilter, mapqFilter,
                                 LOGICAL(reverseComplement)[0],
                                 NA_INTEGER, 0, 0, '\0', '\0', (void *) sbd);
    bd->irange = 0;             /* everything parsed to 'irange' 0 */
    BAM_BUFFER buf = BAMBUFFER(bufext);
    _grow_SCAN_BAM_DATA(bd, buf->n);

    for (int i = 0; i < buf->i; ++i) {
        if (buf->as_mates) {
            sbd->mates_flag = buf->mates[i];
            sbd->partition_id = buf->partition[i];
        }
        int result = _parse1_BAM_DATA(buf->buffer[i], bd);
        if (result < 0) {   /* parse error: e.g., cigar buffer overflow */
            _grow_SCAN_BAM_DATA(bd, 0);
            bd->iparsed = -1;
            break;
        } else if (result == 0L) /* does not pass filter */
            continue;
    }

    int status = bd->iparsed;
    if (status >= 0) {
        _finish1range_BAM_DATA(bd);
        status = bd->iparsed;
    }

    if (status < 0) {
        result = R_NilValue;
        _Free_BAM_DATA(bd);
        UNPROTECT(1);
        Rf_error("bamBuffer 'parse' error code: %d", status);
    }

    _Free_SCAN_BAM_DATA(sbd);
    _Free_BAM_DATA(bd);
    UNPROTECT(1);

    return result;
}

SEXP bambuffer_write(SEXP bufext, SEXP bamext, SEXP filter)
{
    BAM_BUFFER buf;
    BAM_FILE bfile;
    int n, status, filt_n;

    _checkext(bufext, BAMBUFFER_TAG, "bamBuffer 'write'");
    buf = BAMBUFFER(bufext);
    filt_n = Rf_length(filter);
    if ((!IS_LOGICAL(filter)) || !((filt_n == buf->i) || filt_n == 1))
        Rf_error("'filterBam' expected logical(1) or logical(%d)", buf->i);
    _check_isbamfile(bamext, "bamBuffer, 'write'");
    bfile = BAMFILE(bamext);

    n = buf->i;
    for (int i = 0; i < n; ++i)
        if (LOGICAL(filter)[i % filt_n]) {
            status = samwrite(bfile->file, buf->buffer[i]);
            if (status <= 0)
                Rf_error("'bamBuffer' write failed, record %d", i);
        }

    return ScalarInteger(n);
}

SEXP bambuffer_reset(SEXP bufext)
{
    _checkext(bufext, BAMBUFFER_TAG, "bamBuffer 'reset'");
    BAM_BUFFER buf = BAMBUFFER(bufext);
    _bambuffer_reset(buf);
    return ScalarLogical(TRUE);
}
