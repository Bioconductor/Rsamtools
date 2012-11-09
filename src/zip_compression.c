#include <stdio.h>
#include <errno.h>
#include "zip_compression.h"
#include "bgzf.h"
#include "razf.h"

void _zip_error(const char *txt, const char *err, int infd, int outfd)
{
    close(infd);
    close(outfd);
    err ? Rf_error(txt, err) : Rf_error(txt);
}

void _zip_open(SEXP file, SEXP dest, int *infd, int *outfd)
{
    int oflag = O_WRONLY | O_CREAT | O_TRUNC;
#ifdef _WIN32
    oflag |= O_BINARY;
#endif

    if (!IS_CHARACTER(file) || 1L != Rf_length(file))
        Rf_error("'file' must be character(1)");
    if (!IS_CHARACTER(dest) || 1L != Rf_length(dest))
        Rf_error("'dest' must be character(1)");

    *infd = open(translateChar(STRING_ELT(file, 0)), O_RDONLY);
    if (0 > *infd)
        Rf_error("opening 'file': %s", strerror(errno));

    /* we overwrite existing files here */
    *outfd = open(translateChar(STRING_ELT(dest, 0)), oflag, 0666);
    if (0 > *outfd) {
        close(*infd);
        Rf_error("opening 'dest': %s", strerror(errno));
    }
}

void _zip_close(int infd, int outfd)
{
    if (-1L == close(infd))
        _zip_error("closing input after compression: %s",
                   strerror(errno), infd, outfd);
    if (outfd >= 0) {
	if (-1L == close(outfd))
	    Rf_error("closing output after compression: %s",
		     strerror(errno));
    }
}

SEXP bgzip(SEXP file, SEXP dest)
{
    static const int BUF_SIZE = 64 * 1024;
    void *buffer;
    int infd, outfd, cnt;
    BGZF *outp;

    buffer = R_alloc(BUF_SIZE, sizeof(void *));

    _zip_open(file, dest, &infd, &outfd);
    outp = bgzf_fdopen(outfd, "w");
    if (NULL == outp)
        _zip_error("opening output 'dest'", NULL, infd, outfd);

    while (0 < (cnt = read(infd, buffer, BUF_SIZE)))
        if (0 > bgzf_write(outp, buffer, cnt))
            _zip_error("writing compressed output", NULL, infd, outfd);
    if (0 > cnt)
        _zip_error("reading compressed output: %s",
                   strerror(errno), infd, outfd);

    if (0 > bgzf_close(outp))
        Rf_error("closing compressed output");
#ifdef _USE_KNETFILE
    fclose(outp->x.fpw);
#else
    fclose(outp->file);
#endif
    
    _zip_close(infd, -1);

    return dest;
}

SEXP razip(SEXP file, SEXP dest)
{
    static const int WINDOW_SIZE = 4096;
    void *buffer;
    int infd, outfd, cnt;
    RAZF *outp;

    _zip_open(file, dest, &infd, &outfd);

    outp = razf_dopen(outfd, "w");
    if (NULL == outp)
        _zip_error("opening output 'dest'", NULL, infd, outfd);

    buffer = R_alloc(WINDOW_SIZE, sizeof(const int));
    while (0 < (cnt = read(infd, buffer, WINDOW_SIZE)))
        if (0 > razf_write(outp, buffer, cnt))
            _zip_error("writing compressed output", NULL, infd, outfd);
    if (0 > cnt)
        _zip_error("reading compressed output: %s",
                   strerror(errno), infd, outfd);

    razf_close(outp);
    _zip_close(infd, outfd);

    return dest;
}
