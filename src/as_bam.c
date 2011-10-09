#include "samtools/sam.h"
#include "bamfile.h"

int _as_bam(samfile_t * fin, samfile_t * fout)
{
    bam1_t *b = bam_init1();
    int r, count = 0;

    while (0 <= (r = samread(fin, b))) {
        samwrite(fout, b);
        count++;
    }
    bam_destroy1(b);

    return r >= -1 ? count : -1 * count;
}

SEXP as_bam(SEXP file, SEXP destination)
{
    if (!IS_CHARACTER(file) || 1 != LENGTH(file))
        Rf_error("'file' must be character(1)");
    if (!IS_CHARACTER(destination) || 1 != LENGTH(destination))
        Rf_error("'destination' must be character(1)");

    samfile_t *fin =
        _bam_tryopen(translateChar(STRING_ELT(file, 0)), "r", NULL);
    if (fin->header == 0) {
        samclose(fin);
        Rf_error("invalid header");
    }

    samfile_t *fout =           /* f_in leaks if this fails */
        _bam_tryopen(translateChar(STRING_ELT(destination, 0)),
                     "wb", fin->header);

    int count = _as_bam(fin, fout);

    samclose(fin);
    samclose(fout);
    if (count < 0)
        Rf_error("truncated input file at record %d", -1 * count);

    return destination;
}
