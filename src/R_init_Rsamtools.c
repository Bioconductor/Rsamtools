#include <R_ext/Rdynload.h>
#include "io_sam.h"
#include "cigar.h"

#ifdef _WIN32
#include "samtools/knetfile.h"
#endif

static const R_CallMethodDef callMethods[] = {
	/* io_sam.c */
	{".read_bam_header", (DL_FUNC) &read_bam_header, 2},
	{".scan_bam_template", (DL_FUNC) &scan_bam_template, 0},
	{".scan_bam", (DL_FUNC) &scan_bam, 7},
	{".scan_bam_cleanup", (DL_FUNC) &scan_bam_cleanup, 0},
	{".filter_bam", (DL_FUNC) & filter_bam, 8},
	{".count_bam", (DL_FUNC) &count_bam, 6},
	{".index_bam", (DL_FUNC) &index_bam, 1},
	/* cigar.c */
	{".cigar_run_count", (DL_FUNC) &cigar_run_count, 1},
	{".cigar_table", (DL_FUNC) &cigar_table, 1},
	{".split_cigar", (DL_FUNC) &split_cigar, 1},
	{".cigar_to_read_width", (DL_FUNC) &cigar_to_read_width, 2},
	{".cigar_to_IRanges", (DL_FUNC) &cigar_to_IRanges, 3},
	{".cigar_to_list_of_IRanges", (DL_FUNC) &cigar_to_list_of_IRanges, 6},
    {NULL, NULL, 0}
};

void
R_init_Rsamtools(DllInfo *info)
{
#ifdef _WIN32
	int status = knet_win32_init();
	if (status != 0)
	  Rf_error("internal: failed to initialize Winsock; error %d",
			   status);
#endif
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

void
R_unload_Rsamtools(DllInfo *info)
{
#ifdef _WIN32
	knet_win32_destroy();
#endif
}
