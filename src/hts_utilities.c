#include "hts_utilities.h"
#include <cram/cram.h>
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

/*
 * from htslib-1.15.1/sam.c
 *
 * cram_ptell is a pseudo-tell function, because it matches the
 *   position of the disk cursor only after a fresh seek
 *   call. Otherwise it indicates that the read takes place inside the
 *   buffered container previously fetched. It was designed like this
 *   to integrate with the functionality of the iterator stepping
 *   logic.
 */

static int64_t cram_ptell(void *fp)
{
    cram_fd *fd = (cram_fd *)fp;
    cram_container *c;
    cram_slice *s;
    int64_t ret = -1L;

    if (fd) {
        if ((c = fd->ctr) != NULL) {
            if ((s = c->slice) != NULL && s->max_rec) {
                if ((c->curr_slice + s->curr_rec/s->max_rec) >= (c->max_slice + 1))
                    fd->curr_position += c->offset + c->length;
            }
        }
        ret = fd->curr_position;
    }

    return ret;
}

int64_t _hts_utilities_tell(htsFile *fd)
{
    int64_t position = -1;
    if (fd->is_bgzf) {
        position = bgzf_tell(fd->fp.bgzf);
    } else if (fd->is_cram) {
        position = cram_ptell(fd->fp.cram);
    } else {
        Rf_error("[internal] _hts_utilities_tell unknown file format");
    }

    if (position < 0)
        Rf_error("[internal] _hts_utilities_tell failed");

    return position;
}
