#include <R_ext/Rdynload.h>
#include "io_sam.h"
#include "cigar.h"
#include "Alignments_class.h"

#ifdef _WIN32
#include "samtools/knetfile.h"
#endif

static const R_CallMethodDef callMethods[] = {
	/* io_sam.c */
	{".read_bam_header", (DL_FUNC) &read_bam_header, 2},
	{".scan_bam_template", (DL_FUNC) &scan_bam_template, 0},
	{".scan_bam", (DL_FUNC) &scan_bam, 8},
	{".scan_bam_cleanup", (DL_FUNC) &scan_bam_cleanup, 0},
	{".filter_bam", (DL_FUNC) & filter_bam, 8},
	{".count_bam", (DL_FUNC) &count_bam, 6},
	{".index_bam", (DL_FUNC) &index_bam, 1},
	/* cigar.c */
	{".cigar_run_count", (DL_FUNC) &cigar_run_count, 1},
	{".cigar_table", (DL_FUNC) &cigar_table, 1},
	{".valid_cigar", (DL_FUNC) &valid_cigar, 2},
	{".split_cigar", (DL_FUNC) &split_cigar, 1},
	{".cigar_op_table", (DL_FUNC) &cigar_op_table, 1},
	{".cigar_to_qwidth", (DL_FUNC) &cigar_to_qwidth, 2},
	{".cigar_to_qwidth2", (DL_FUNC) &cigar_to_qwidth2, 2},
	{".cigar_to_width", (DL_FUNC) &cigar_to_width, 1},
	{".cigar_qnarrow", (DL_FUNC) &cigar_qnarrow, 3},
	{".cigar_to_IRanges", (DL_FUNC) &cigar_to_IRanges, 3},
	{".cigar_to_list_of_IRanges_by_alignment",
	 (DL_FUNC) &cigar_to_list_of_IRanges_by_alignment, 4},
	{".cigar_to_list_of_IRanges_by_rname",
	 (DL_FUNC) &cigar_to_list_of_IRanges_by_rname, 6},
	/* Alignments_class.c */
	{".logical_as_compact_raw_vector", (DL_FUNC) &logical_as_compact_raw_vector, 1},
	{".compact_raw_vector_as_logical", (DL_FUNC) &compact_raw_vector_as_logical, 2},
	{".subset_compact_raw_vector", (DL_FUNC) &subset_compact_raw_vector, 2},
	{".compact_raw_vector_bit_count", (DL_FUNC) &compact_raw_vector_bit_count, 1},
	{".compact_raw_vector_last_bit", (DL_FUNC) &compact_raw_vector_last_bit, 1},
	{".compact_raw_vector_set_op", (DL_FUNC) &compact_raw_vector_set_op, 3},
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
