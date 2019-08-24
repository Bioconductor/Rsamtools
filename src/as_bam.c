#include <samtools-1.7-compat.h>
#include "bamfile.h"

int _as_bam(htsFile * fin, htsFile * fout, bam_hdr_t *header)
{
    bam1_t *b = bam_init1();
    int r, count = 0;

    while (0 <= (r = sam_read1(fin, header, b))) {
        if (sam_write1(fout, header, b) == 0)
            break;
        count++;
    }
    bam_destroy1(b);

    return r >= -1 ? count : -1 * count;
}

SEXP as_bam(SEXP file, SEXP destination, SEXP binary)
{
    if (!IS_CHARACTER(file) || 1 != LENGTH(file))
        Rf_error("'file' must be character(1)");
    if (!IS_CHARACTER(destination) || 1 != LENGTH(destination))
        Rf_error("'destination' must be character(1)");
    if (!IS_LOGICAL(binary) || 1 != LENGTH(binary))
        Rf_error("'binary' must be logical(1)");

    const char
        *file_c = translateChar(STRING_ELT(file, 0)),
        *destination_c = translateChar(STRING_ELT(destination, 0));
    const char *mode = LOGICAL(binary)[0] ? "wb" : "w";

    htsFile *fin, *fout;
    bam_hdr_t *header;

    fin = _bam_tryopen(file_c, "r", NULL);
    header = sam_hdr_read(fin);
    if (header == 0) {
        sam_close(fin);
        Rf_error("invalid header");
    }

    fout = _bam_tryopen(destination_c, mode, NULL);
    if (sam_hdr_write(fout, header) < 0) {
        bam_hdr_destroy(header);
        sam_close(fin);
        sam_close(fout);
        Rf_error("failed to write header to output file");
    }

    int count = _as_bam(fin, fout, header);

    if (sam_close(fin) < 0)
        Rf_warning("'asBam' failed to close input file");
    if (sam_close(fout) < 0)
        Rf_warning("'asBam' failed to closue output file");
    if (count < 0)
        Rf_error("truncated input file at record %d", -1 * count);

    return destination;
}
