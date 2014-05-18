#include "bamfile.h"
#include "io_sam.h"
#include "bam_mate_iter.h"
#include "utilities.h"

static SEXP BAMFILE_TAG = NULL;
#define TYPE_BAM 1

void _check_isbamfile(SEXP ext, const char *lbl)
{
    _checkext(ext, BAMFILE_TAG, lbl);
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
    BAM_FILE bfile = BAMFILE(ext);
    if (NULL != bfile->file)
        samclose(bfile->file);
    if (NULL != bfile->index)
        bam_index_destroy(bfile->index);
    if (NULL != bfile->iter)
        bam_mate_iter_destroy(bfile->iter);
    bfile->file = NULL;
    bfile->index = NULL;
    bfile->iter = NULL;
}

static void _bamfile_finalizer(SEXP ext)
{
    if (NULL == R_ExternalPtrAddr(ext))
        return;
    _bamfile_close(ext);
    BAM_FILE bfile = BAMFILE(ext);
    Free(bfile);
    R_SetExternalPtrAddr(ext, NULL);
}

SEXP bamfile_init()
{
    BAMFILE_TAG = install("BamFile");
    return R_NilValue;
}

static BAM_FILE _bamfile_open_r(SEXP filename, SEXP indexname, SEXP filemode)
{
    BAM_FILE bfile = (BAM_FILE) Calloc(1, _BAM_FILE);

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
        bfile->irange0 = 0;
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

    bfile->iter = NULL;
    return bfile;
}

static BAM_FILE _bamfile_open_w(SEXP file0, SEXP file1)
{
    samfile_t *infile, *outfile;
    BAM_FILE bfile;

    if (0 == Rf_length(file1))
        Rf_error("'file1' must be a character(1) path to a valid bam file");
    infile = _bam_tryopen(translateChar(STRING_ELT(file1, 0)), "rb", 0);
    outfile = _bam_tryopen(translateChar(STRING_ELT(file0, 0)), "wb",
                           infile->header);
    samclose(infile);

    bfile = (BAM_FILE) Calloc(1, _BAM_FILE);
    bfile->file = outfile;
    bfile->pos0 = bam_tell(bfile->file->x.bam);
    bfile->irange0 = 0;

    return bfile;
}

SEXP bamfile_open(SEXP file0, SEXP file1, SEXP mode)
{
    _checknames(file0, file1, mode);
    BAM_FILE bfile;
    if (*CHAR(STRING_ELT(mode, 0)) == 'r')
        bfile = _bamfile_open_r(file0, file1, mode);
    else
        bfile = _bamfile_open_w(file0, file1);

    SEXP ext = PROTECT(R_MakeExternalPtr(bfile, BAMFILE_TAG, file0));
    R_RegisterCFinalizerEx(ext, _bamfile_finalizer, TRUE);
    UNPROTECT(1);

    return ext;
}

SEXP bamfile_close(SEXP ext)
{
    _checkext(ext, BAMFILE_TAG, "close");
    _bamfile_close(ext);
    return ext;
}

SEXP bamfile_isopen(SEXP ext)
{
    int ans = FALSE;
    if (NULL != BAMFILE(ext)) {
        _checkext(ext, BAMFILE_TAG, "isOpen");
        ans = NULL != BAMFILE(ext)->file;
    }
    return ScalarLogical(ans);
}

SEXP bamfile_isincomplete(SEXP ext)
{
    int ans = FALSE;
    BAM_FILE bfile;
    if (NULL != BAMFILE(ext)) {
        _checkext(ext, BAMFILE_TAG, "isIncomplete");
        bfile = BAMFILE(ext);
        if (NULL != bfile && NULL != bfile->file) {
            /* heuristic: can we read a record? bam_seek does not
             * support SEEK_END */
            off_t offset = bam_tell(bfile->file->x.bam);
            char buf;
            ans = bam_read(bfile->file->x.bam, &buf, 1) > 0;
            bam_seek(bfile->file->x.bam, offset, SEEK_SET);
        }
    }
    return ScalarLogical(ans);
}

/* implementation */

SEXP read_bamfile_header(SEXP ext, SEXP what)
{
    _checkext(ext, BAMFILE_TAG, "scanBamHeader");
    if (!(IS_LOGICAL(what) && (2L == LENGTH(what))))
        Rf_error("'what' must be logical(2)");
    if (!LOGICAL(bamfile_isopen(ext))[0])
        Rf_error("open() BamFile before reading header");
    return _read_bam_header(ext, what);
}

SEXP scan_bamfile(SEXP ext, SEXP space, SEXP keepFlags, SEXP isSimpleCigar,
                  SEXP reverseComplement, SEXP yieldSize, SEXP template_list,
                  SEXP obeyQname, SEXP asMates)
{
    _checkext(ext, BAMFILE_TAG, "scanBam");
    _checkparams(space, keepFlags, isSimpleCigar);
    if (!(IS_LOGICAL(reverseComplement) && (1L == LENGTH(reverseComplement))))
        Rf_error("'reverseComplement' must be logical(1)");
    if (!(IS_INTEGER(yieldSize) && (1L == LENGTH(yieldSize))))
        Rf_error("'yieldSize' must be integer(1)");
    if (!(IS_LOGICAL(obeyQname) && (1L == LENGTH(obeyQname))))
        Rf_error("'obeyQname' must be logical(1)");
    if (!(IS_LOGICAL(asMates) && (1L == LENGTH(asMates))))
        Rf_error("'asMates' must be logical(1)");
    _bam_check_template_list(template_list);
    return _scan_bam(ext, space, keepFlags, isSimpleCigar,
                     reverseComplement, yieldSize, template_list, 
                     obeyQname, asMates);
}

SEXP count_bamfile(SEXP ext, SEXP space, SEXP keepFlags, SEXP isSimpleCigar)
{
    _checkext(ext, BAMFILE_TAG, "countBam");
    _checkparams(space, keepFlags, isSimpleCigar);
    SEXP count = _count_bam(ext, space, keepFlags, isSimpleCigar);
    if (R_NilValue == count)
        Rf_error("'countBam' failed");
    return count;
}

SEXP prefilter_bamfile(SEXP ext, SEXP space, SEXP keepFlags, SEXP isSimpleCigar,
                       SEXP yieldSize, SEXP obeyQname, SEXP asMates)
{
    _checkext(ext, BAMFILE_TAG, "filterBam");
    _checkparams(space, keepFlags, isSimpleCigar);
    if (!(IS_INTEGER(yieldSize) && (1L == LENGTH(yieldSize))))
        Rf_error("'yieldSize' must be integer(1)");
    if (!(IS_LOGICAL(obeyQname) && (1L == LENGTH(obeyQname))))
        Rf_error("'obeyQname' must be logical(1)");
    if (!(IS_LOGICAL(asMates) && (1L == LENGTH(asMates))))
        Rf_error("'asMates' must be logical(1)");
    SEXP result =
        _prefilter_bam(ext, space, keepFlags, isSimpleCigar, yieldSize,
                       obeyQname, asMates);
    if (R_NilValue == result)
        Rf_error("'filterBam' failed during pre-filtering");
    return result;
}

SEXP filter_bamfile(SEXP ext, SEXP space, SEXP keepFlags, SEXP isSimpleCigar,
                    SEXP fout_name, SEXP fout_mode)
{
    _checkext(ext, BAMFILE_TAG, "filterBam");
    _checkparams(space, keepFlags, isSimpleCigar);
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

