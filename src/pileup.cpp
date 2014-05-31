#include "io_sam.h"
#include "utilities.h"
#include "PileupBufferShim.h"
#include "pileup.h"

static int _filter_and_parse1_pileup(const bam1_t *bam, void *data)
{
    BAM_DATA bd = (BAM_DATA) data;
    int result = _filter1_BAM_DATA(bam, bd);
    if (result < 0) {
        // possible parse error: e.g., cigar buf overflow
        Rf_error("parsing BAM file failed. Is the file corrupt?");
    }
    else if(result) {
        PileupBufferShim *shim = (PileupBufferShim *) bd->extra;
        // push to buffer the alignments that meet criteria
        shim->plbuf_push(bam);
    }
    bd->iparsed += 1;
    return result;
}

static void _finish1range_pileup(BAM_DATA bd)
{
    // invoked each time a range is finished
    PileupBufferShim *shim = (PileupBufferShim *) bd->extra;
    shim->finish1(bd->irange);
    if (bd->irange + 1 < bd->nrange)
        shim->start1(bd->irange + 1);
}

static SEXP _pileup_bam_result_init(SEXP space)
{
    // space is list of vectors, so length of outer list
    // is num ranges
    const int nrange =
        R_NilValue == space ? 1 : Rf_length(VECTOR_ELT(space, 0));
    return Rf_allocVector(VECSXP, nrange);
}

static SEXP _pileup_bam(SEXP ext, SEXP space, SEXP keepFlags,
    SEXP reverseComplement, SEXP isSimpleCigar, SEXP yieldSize,
    SEXP obeyQname, SEXP asMates, SEXP qnamePrefixEnd,
    SEXP qnameSuffixStart, PileupBuffer& buffer)
{
    _check_isbamfile(ext, "pileup");
    _checkparams(space, keepFlags, isSimpleCigar);
    if (!(IS_INTEGER(yieldSize) && (1L == LENGTH(yieldSize))))
        Rf_error("'yieldSize' must be integer(1)");
    if (!(IS_LOGICAL(obeyQname) && (1L == LENGTH(obeyQname))))
        Rf_error("'obeyQname' must be logical(1)");
    if (!(IS_LOGICAL(asMates) && (1L == LENGTH(asMates))))
        Rf_error("'asMates' must be logical(1)");

    SEXP result = PROTECT(_pileup_bam_result_init(space));
    PileupBufferShim shim(space, result, buffer);
    char qname_prefix = '\0';
    SEXP prefix_elt = STRING_ELT(qnamePrefixEnd, 0);
    if (prefix_elt != NA_STRING);
        qname_prefix = CHAR(prefix_elt)[0];
    char qname_suffix = '\0';
    SEXP suffix_elt = STRING_ELT(qnameSuffixStart, 0);
    if (suffix_elt != NA_STRING);
        qname_suffix = CHAR(suffix_elt)[0];
    BAM_DATA bd = _init_BAM_DATA(ext, space, keepFlags, isSimpleCigar,
                                 LOGICAL(reverseComplement)[0],
                                 INTEGER(yieldSize)[0],
                                 LOGICAL(obeyQname)[0], 
                                 LOGICAL(asMates)[0],
                                 qname_prefix, qname_suffix,
                                 (void *) &shim);
    int status = 0;
    if (bd->irange < bd->nrange) {
        shim.start1(bd->irange);
        // _do_scan_bam found in io_sam.c;
        // mostly a wrapper for
        // _scan_bam_fetch(BAM_DATA bd, SEXP space, int* start, int* end,
        //                 bam_fetch_f , bam_fetch_mate_f , _FINISH_FUNC)
        status = _do_scan_bam(bd, space, _filter_and_parse1_pileup,
                              NULL, _finish1range_pileup);
    }
    
    if (status < 0) {
        int idx = bd->irec;
        int parse_status = bd->parse_status;
        _Free_BAM_DATA(bd);
        status = parse_status;
        Rf_error("'pileup' failed:\n  record: %d\n  error: %d",
                 idx, parse_status);
    }

    _Free_BAM_DATA(bd);
    UNPROTECT(1);
    return result;
}

extern "C" {
    SEXP c_Pileup(SEXP ext, SEXP space, SEXP keepFlags,
                  SEXP isSimpleCigar, SEXP reverseComplement,
                  SEXP yieldSize, SEXP obeyQname, SEXP asMates,
                  SEXP qnamePrefixEnd, SEXP qnameSuffixStart, 
                  SEXP schema, SEXP pileupParams, SEXP seqnamesLevels)
    {
        if (!IS_LIST(schema))
            Rf_error("'schema' must be list()");
        if (!IS_LIST(pileupParams))
            Rf_error("'pileupParams' must be list()");

        Pileup buffer = Pileup(schema, pileupParams, seqnamesLevels);
        return _pileup_bam(ext, space, keepFlags,
            reverseComplement, isSimpleCigar, yieldSize, obeyQname,
            asMates, qnamePrefixEnd, qnameSuffixStart, buffer);
    }
}
