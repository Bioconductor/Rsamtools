#include "hts_utilities.h"
#include <htslib/cram.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <R_ext/Error.h>

htsFormat _hts_utilities_format(const char *filename)
{
    hFILE *fp = hopen(filename, "r");
    if (fp == NULL)
        Rf_error("can't open file");
    htsFormat fmt;
    if (hts_detect_format(fp, &fmt) < 0) {
        hclose_abruptly(fp);
        Rf_error("failed to detect file format");
    }
    if (hclose(fp) < 0)
        Rf_error("failed to close file after detecting format");

    return fmt;
}

int _hts_utilities_seek(htsFile *fd, off_t offset, int whence)
{
    int status = -1;

    if (fd->is_bgzf) {
        status = bgzf_seek(fd->fp.bgzf, offset, whence);
    } else if (fd->is_cram) {
        status = cram_seek(fd->fp.cram, offset, whence);
    } else {
        Rf_error("[internal] _hts_utilities_seek unknown file format");
    }

    if (status < 0)
        Rf_error("[internal] _hts_utilities_seek failed");

    return status;
}

int64_t _hts_utilities_tell(htsFile *fd)
{
    int64_t position = -1;
    if (fd->is_bgzf) {
        position = bgzf_tell(fd->fp.bgzf);
    } else if (fd->is_cram) {
        position = cram_tell(fd->fp.cram);
    } else {
        Rf_error("[internal] _hts_utilities_tell unknown file format");
    }

    if (position < 0)
        Rf_error("[internal] _hts_utilities_tell failed");

    return position;
}
