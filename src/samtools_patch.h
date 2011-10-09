#ifndef SAMTOOLS_PATCH_H
#define SAMTOOLS_PATCH_H

/* capture samtools errors */
extern void _samtools_abort();
extern void _samtools_exit(int status);
extern int _samtools_fprintf(FILE *, const char *, ...);

#endif
