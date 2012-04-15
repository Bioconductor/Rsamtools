#include "bamfile.h"
#include "io_sam.h"
#include "utilities.h"

static SEXP BAMFILE_TAG = NULL;
#define TYPE_BAM 1

void _check_isbamfile(SEXP ext, const char *lbl)
{
    _scan_checkext(ext, BAMFILE_TAG, lbl);
}

samfile_t *_bam_tryopen(const char *filename, const char *filemode, void *aux)
{
    samfile_t *sfile = samopen(filename, filemode, aux);
    if (sfile == 0)
        Rf_error("failed to open SAM/BAM file\n  file: '%s'", filename);
    if (sfile->header == 0 || sfile->header->n_targets == 0) {
        samclose(sfile);
        Rf_error("SAM/BAM header missing or empty\n  file: '%s'", filename);
    }
    return sfile;
}

static bam_index_t *_bam_tryindexload(const char *indexname)
{
    bam_index_t *index = bam_index_load(indexname);
    if (index == 0)
        Rf_error("failed to load BAM index\n  file: %s", indexname);
    return index;
}

static void _bamfile_close(SEXP ext)
{
    _BAM_FILE *bfile = BAMFILE(ext);
    if (NULL != bfile->file)
        samclose(bfile->file);
    if (NULL != bfile->index)
        bam_index_destroy(bfile->index);
    bfile->file = NULL;
    bfile->index = NULL;
}

static void _bamfile_finalizer(SEXP ext)
{
    if (NULL == R_ExternalPtrAddr(ext))
        return;
    _bamfile_close(ext);
    _BAM_FILE *bfile = BAMFILE(ext);
    Free(bfile);
    R_SetExternalPtrAddr(ext, NULL);
}

SEXP bamfile_init()
{
    BAMFILE_TAG = install("BamFile");
    return R_NilValue;
}

SEXP bamfile_open(SEXP filename, SEXP indexname, SEXP filemode)
{
    _scan_checknames(filename, indexname, filemode);

    _BAM_FILE *bfile = (_BAM_FILE *) Calloc(1, _BAM_FILE);

    bfile->file = NULL;
    if (0 != Rf_length(filename)) {
        const char *cfile = translateChar(STRING_ELT(filename, 0));
        bfile->file = _bam_tryopen(cfile, CHAR(STRING_ELT(filemode, 0)), 0);
        if ((bfile->file->type & TYPE_BAM) != 1) {
            samclose(bfile->file);
            Free(bfile);
            Rf_error("'filename' is not a BAM file\n  file: %s", cfile);
        }
        bfile->pos0 = bam_tell(bfile->file->x.bam);
    }

    bfile->index = NULL;
    if (0 != Rf_length(indexname)) {
        const char *cindex = translateChar(STRING_ELT(indexname, 0));
        bfile->index = _bam_tryindexload(cindex);
        if (NULL == bfile->index) {
            samclose(bfile->file);
            Free(bfile);
            Rf_error("failed to open BAM index\n  index: %s\n", cindex);
        }
    }

    SEXP ext = PROTECT(R_MakeExternalPtr(bfile, BAMFILE_TAG, filename));
    R_RegisterCFinalizerEx(ext, _bamfile_finalizer, TRUE);
    UNPROTECT(1);

    return ext;
}

SEXP bamfile_close(SEXP ext)
{
    _scan_checkext(ext, BAMFILE_TAG, "close");
    _bamfile_close(ext);
    return ext;
}

SEXP bamfile_isopen(SEXP ext)
{
    SEXP ans = ScalarLogical(FALSE);
    if (NULL != BAMFILE(ext)) {
        _scan_checkext(ext, BAMFILE_TAG, "isOpen");
        if (NULL != BAMFILE(ext)->file)
            ans = ScalarLogical(TRUE);
    }
    return ans;
}

/* implementation */

SEXP read_bamfile_header(SEXP ext)
{
    _scan_checkext(ext, BAMFILE_TAG, "scanBamHeader");
    return _read_bam_header(ext);
}

SEXP scan_bamfile(SEXP ext, SEXP space, SEXP keepFlags, SEXP isSimpleCigar,
                  SEXP reverseComplement, SEXP yieldSize, SEXP template_list)
{
    _scan_checkext(ext, BAMFILE_TAG, "scanBam");
    _scan_checkparams(space, keepFlags, isSimpleCigar);
    if (!(IS_LOGICAL(reverseComplement) && (1L == LENGTH(reverseComplement))))
        Rf_error("'reverseComplement' must be logical(1)");
    if (!(IS_INTEGER(yieldSize) && (1L == LENGTH(yieldSize))))
        Rf_error("'yieldSize' must be integer(1)");
    _bam_check_template_list(template_list);
    return _scan_bam(ext, space, keepFlags, isSimpleCigar,
                     reverseComplement, yieldSize, template_list);
}

SEXP count_bamfile(SEXP ext, SEXP space, SEXP keepFlags, SEXP isSimpleCigar)
{
    _scan_checkext(ext, BAMFILE_TAG, "countBam");
    _scan_checkparams(space, keepFlags, isSimpleCigar);
    SEXP count = _count_bam(ext, space, keepFlags, isSimpleCigar);
    if (R_NilValue == count)
        Rf_error("'countBam' failed");
    return count;
}

SEXP filter_bamfile(SEXP ext, SEXP space, SEXP keepFlags, SEXP isSimpleCigar,
                    SEXP fout_name, SEXP fout_mode)
{
    _scan_checkext(ext, BAMFILE_TAG, "filterBam");
    _scan_checkparams(space, keepFlags, isSimpleCigar);
    if (!IS_CHARACTER(fout_name) || 1 != LENGTH(fout_name))
        Rf_error("'fout_name' must be character(1)");
    if (!IS_CHARACTER(fout_mode) || 1 != LENGTH(fout_mode))
        Rf_error("'fout_mode' must be character(1)");
    SEXP result = _filter_bam(ext, space, keepFlags, isSimpleCigar,
                              fout_name, fout_mode);
    if (R_NilValue == result)
        Rf_error("'filterBam' failed");
    return result;
}
