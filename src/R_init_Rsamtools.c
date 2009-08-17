#include "Rsamtools.h"

static const R_CallMethodDef callMethods[] = {
	/* io_sam.c */
	{".read_bam_header", (DL_FUNC) &read_bam_header, 3},
	{".scan_bam_open", (DL_FUNC) &scan_bam_open, 2},
	{".scan_bam_template", (DL_FUNC) &scan_bam_template, 0},
	{".scan_bam", (DL_FUNC) &scan_bam, 5},
    {NULL, NULL, 0}
};

void
R_init_Rsamtools(DllInfo *info)
{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

