#include <stdio.h>
#include <errno.h>
#include <zlib.h>
#include <fcntl.h>
#include <htslib/bgzf.h>
#include "zip_compression.h"

static void _zip_error(const char *txt, const char *err, int infd, int outfd)
{
    close(infd);
    close(outfd);
    err ? Rf_error(txt, err) : Rf_error("%s", txt);
}

static void _zip_open(SEXP file, SEXP dest, int *infd, int *outfd)
{
    int iflag = O_RDONLY, oflag = O_WRONLY | O_CREAT | O_TRUNC;
#ifdef _WIN32
    iflag |= O_BINARY;
    oflag |= O_BINARY;
#endif

    if (!IS_CHARACTER(file) || 1L != Rf_length(file))
        Rf_error("'file' must be character(1)");
    if (!IS_CHARACTER(dest) || 1L != Rf_length(dest))
        Rf_error("'dest' must be character(1)");

    *infd = open(translateChar(STRING_ELT(file, 0)), iflag);
    if (0 > *infd)
        Rf_error("opening 'file': %s", strerror(errno));

    /* we overwrite existing files here */
    *outfd = open(translateChar(STRING_ELT(dest, 0)), oflag, 0666);
    if (0 > *outfd) {
        close(*infd);
        Rf_error("opening 'dest': %s", strerror(errno));
    }
}

SEXP bgzip(SEXP file, SEXP dest)
{
    static const int BUF_SIZE = 64 * 1024;
    void *buffer;
    int infd, outfd, cnt;
    gzFile in;
    BGZF *outp;

    buffer = R_alloc(BUF_SIZE, sizeof(void *));

    _zip_open(file, dest, &infd, &outfd);
    in = gzdopen(infd, "rb");
    if (NULL == in)
        _zip_error("opening input 'file'", NULL, infd, outfd);
    outp = bgzf_dopen(outfd, "w");
    if (NULL == outp)
        _zip_error("opening output 'dest'", NULL, infd, outfd);

    while (0 < (cnt = gzread(in, buffer, BUF_SIZE)))
        if (0 > bgzf_write(outp, buffer, cnt))
            _zip_error("writing compressed output", NULL, infd, outfd);
    if (0 > cnt)
        _zip_error("reading compressed input: %s",
                   strerror(errno), infd, outfd);

    if (0 > bgzf_close(outp))
        Rf_error("closing compressed output");
    if (gzclose(in) != Z_OK)
        _zip_error("closing input after compression", NULL, infd, outfd);

    return dest;
}

