#include <R_ext/Rdynload.h>
#include "zip_compression.h"
#include "bamfile.h"
#include "bcffile.h"
#include "fafile.h"
#include "tabixfile.h"
#include "io_sam.h"
#include "as_bam.h"
#include "pileupbam.h"

#ifdef _WIN32
#include "samtools/knetfile.h"
#endif

static const R_CallMethodDef callMethods[] = {
    /* zip_compression.c */
    {".bgzip", (DL_FUNC) & bgzip, 2},
    {".razip", (DL_FUNC) & razip, 2},
    /* bamfile.c */
    {".bamfile_init", (DL_FUNC) & bamfile_init, 0},
    {".bamfile_open", (DL_FUNC) & bamfile_open, 3},
    {".bamfile_close", (DL_FUNC) & bamfile_close, 1},
    {".bamfile_isopen", (DL_FUNC) & bamfile_isopen, 1},
    {".read_bamfile_header", (DL_FUNC) & read_bamfile_header, 1},
    {".scan_bamfile", (DL_FUNC) & scan_bamfile, 7},
    {".count_bamfile", (DL_FUNC) & count_bamfile, 4},
    {".filter_bamfile", (DL_FUNC) & filter_bamfile, 6},
    /* as_bam.c */
    {".as_bam", (DL_FUNC) & as_bam, 2},
    /* io_sam.c */
    {".scan_bam_template", (DL_FUNC) & scan_bam_template, 1},
    {".scan_bam_cleanup", (DL_FUNC) & scan_bam_cleanup, 0},
    {".sort_bam", (DL_FUNC) & sort_bam, 4},
    {".merge_bam", (DL_FUNC) & merge_bam, 8},
    {".index_bam", (DL_FUNC) & index_bam, 1},
    /* bcffile.c */
    {".bcffile_init", (DL_FUNC) & bcffile_init, 0},
    {".bcffile_open", (DL_FUNC) & bcffile_open, 3},
    {".bcffile_close", (DL_FUNC) & bcffile_close, 1},
    {".bcffile_isopen", (DL_FUNC) & bcffile_isopen, 1},
    {".bcffile_isvcf", (DL_FUNC) & bcffile_isvcf, 1},
    {".scan_bcf_header", (DL_FUNC) & scan_bcf_header, 1},
    {".scan_bcf", (DL_FUNC) & scan_bcf, 3},
    {".as_bcf", (DL_FUNC) & as_bcf, 3},
    {".index_bcf", (DL_FUNC) & index_bcf, 1},
    /* fafile.c */
    {".fafile_init", (DL_FUNC) & fafile_init, 0},
    {".fafile_open", (DL_FUNC) & fafile_open, 1},
    {".fafile_close", (DL_FUNC) & fafile_close, 1},
    {".fafile_isopen", (DL_FUNC) & fafile_isopen, 1},
    {".index_fa", (DL_FUNC) & index_fa, 1},
    {".n_fa", (DL_FUNC) & n_fa, 1},
    {".scan_fa", (DL_FUNC) & scan_fa, 5},
    /* tabixfile */
    {".tabixfile_init", (DL_FUNC) & tabixfile_init, 0},
    {".tabixfile_open", (DL_FUNC) & tabixfile_open, 2},
    {".tabixfile_close", (DL_FUNC) & tabixfile_close, 1},
    {".tabixfile_isopen", (DL_FUNC) & tabixfile_isopen, 1},
    {".index_tabix", (DL_FUNC) & index_tabix, 8},
    {".header_tabix", (DL_FUNC) & header_tabix, 1},
    {".scan_tabix", (DL_FUNC) & scan_tabix, 5},
    {".tabix_as_character", (DL_FUNC) & tabix_as_character, 4},
    /* pileup */
    {".apply_pileups", (DL_FUNC) & apply_pileups, 5},
    {NULL, NULL, 0}
};

void R_init_Rsamtools(DllInfo * info)
{
#ifdef _WIN32
    int status = knet_win32_init();
    if (status != 0)
        Rf_error("internal: failed to initialize Winsock; error %d", status);
#endif
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

void R_unload_Rsamtools(DllInfo *info)
{
    (void) info;
#ifdef _WIN32
    knet_win32_destroy();
#endif
}
