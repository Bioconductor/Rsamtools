#include <Rdefines.h>
#include "samtools_patch.h"

void _samtools_exit(int status)
{
    Rf_error
        ("internal: samtools invoked 'exit(%d)'; see warnings() and restart R",
         status);
}

void _samtools_abort()
{
    Rf_error
        ("internal: samtools invoked 'abort'; see warnings() and restart R");
}

int _samtools_fprintf(FILE * file, const char *fmt, ...)
{
    static const int bufsize = 2048;
    va_list argp;
    int n;

    if (stderr != file) {
        va_start(argp, fmt);
        n = vfprintf(file, fmt, argp);
        va_end(argp);
    } else {
        /* silence some messages */
        char *buf = (char *) R_alloc(bufsize, sizeof(char));
        if (0 == strncmp("[samopen] SAM header is present:", fmt, 32) ||
            0 == strncmp("[fai_load] build FASTA index.", fmt, 29))
            return 0;
        va_start(argp, fmt);
        n = vsnprintf(buf, bufsize, fmt, argp);
        va_end(argp);
        Rf_warning(buf);
    }
    return n;
}
