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

static SEXP _pileup_bam_result_init(SEXP regions)
{
    // regions is list of vectors, so length of outer list
    // is num ranges
    const int nrange =
        R_NilValue == regions ? 1 : Rf_length(VECTOR_ELT(regions, 0));
    return Rf_allocVector(VECSXP, nrange);
}

static void _finish_buffered_yieldSize_chunk(BAM_DATA bd) {
    PileupBufferShim *shim = static_cast<PileupBufferShim*>(bd->extra);
    shim->process_yieldSize_chunk(); // run Pileup::insert loop
}

static SEXP _pileup_bam(SEXP ext, SEXP regions, SEXP keepFlags,
    SEXP reverseComplement, SEXP isSimpleCigar,
    SEXP tagFilter, SEXP mapqFilter,
    SEXP yieldSize, SEXP obeyQname, SEXP asMates, SEXP qnamePrefixEnd,
    SEXP qnameSuffixStart, PileupBuffer& buffer)
{
    _check_isbamfile(ext, "pileup");
    _checkparams(regions, keepFlags, isSimpleCigar);
    if (!(Rf_isInteger(yieldSize) && (1L == Rf_length(yieldSize))))
        Rf_error("'yieldSize' must be integer(1)");
    if (!(Rf_isLogical(obeyQname) && (1L == Rf_length(obeyQname))))
        Rf_error("'obeyQname' must be logical(1)");
    if (!(Rf_isLogical(asMates) && (1L == Rf_length(asMates))))
        Rf_error("'asMates' must be logical(1)");

    SEXP result = PROTECT(_pileup_bam_result_init(regions));
    PileupBufferShim shim(regions, result, buffer);
    BAM_DATA bd = _init_BAM_DATA(ext, regions, keepFlags, isSimpleCigar,
                                 tagFilter, mapqFilter,
                                 LOGICAL(reverseComplement)[0],
                                 INTEGER(yieldSize)[0],
                                 LOGICAL(obeyQname)[0], 
                                 LOGICAL(asMates)[0], '\0', '\0',
                                 (void *) &shim);
    int status = 0;
    if(!(dynamic_cast<Pileup&>(buffer).isBufferedPileup())) {
        if (bd->irange < bd->nrange) {
            shim.start1(bd->irange);
            // _do_scan_bam found in io_sam.c;
            // mostly a wrapper for
            // _scan_bam_fetch(BAM_DATA bd, SEXP space, int* start, int* end,
            //                 bam_fetch_f , bam_fetch_mate_f , _FINISH_FUNC)
            status = _do_scan_bam(bd, regions, _filter_and_parse1_pileup,
                                  (bam_fetch_mate_f) NULL,
                                  _finish1range_pileup);
        }
    } else { // must do loop if buffered
        shim.start1(0);
        status = _do_scan_bam(bd, regions, _filter_and_parse1_pileup,
                              (bam_fetch_mate_f) NULL,
                              _finish_buffered_yieldSize_chunk);
        for(;;) {
            if(!(dynamic_cast<Pileup&>(buffer).needMoreInput()) || status <= 0)
                break;
            status = _do_scan_bam(bd, regions, _filter_and_parse1_pileup,
                                  (bam_fetch_mate_f) NULL,
                                  _finish_buffered_yieldSize_chunk);
        }
        shim.finish1(0);
    }

    // the decision to deallocate must happen in here before return
    // the SEXP because Rf_error would lead to the dynamic memory
    // being left behind
    if(status <= 0) { // error or EOI; must deallocate dynamic memory
        // rename 'deallocate and flush'? two separate functions (so
        // that deallocate can be called without flushing?)
        buffer.signalEOI();
        shim.flush();
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

static SEXP _bamheaderAsSeqnames(bam_header_t *header) {
    if(header == NULL)
        Rf_error("'header' must not be NULL");
    int numSeqnames = header->n_targets;
    SEXP seqnames = PROTECT(Rf_allocVector(STRSXP, numSeqnames));
    for(int i = 0; i != numSeqnames; ++i) {
        SET_STRING_ELT(seqnames, i, Rf_mkChar(header->target_name[i]));
    }

    UNPROTECT(1);
    return seqnames;
}

class PosCacheColl;

extern "C" {
    SEXP c_Pileup(SEXP ext, SEXP regions, SEXP keepFlags,
                  SEXP isSimpleCigar, SEXP tagFilter, SEXP mapqFilter,
                  SEXP reverseComplement, SEXP yieldSize,
                  SEXP obeyQname, SEXP asMates,
                  SEXP qnamePrefixEnd, SEXP qnameSuffixStart, 
                  SEXP pileupParams)
    {
        if (!Rf_isVector(pileupParams))
            Rf_error("'pileupParams' must be list()");
        SEXP seqnamesLevels =
            PROTECT(_bamheaderAsSeqnames(BAMFILE(ext)->file->header));
        // 'ranged' means user asked for specific genomic range(s) by
        // providing 'which' argument to 'ScanBamParam'
        bool isRanged = regions != R_NilValue;

        // This check has to happen at this level; can't use value of
        // yieldSize further down in hierarchy
        // 'buffered' means pileup will store pileup results for
        // partially completed genomic positions across calls to
        // pileup
        bool isBuffered = !isRanged && INTEGER(yieldSize)[0] != NA_INTEGER;

        Pileup buffer =
            Pileup(isRanged, isBuffered, pileupParams, seqnamesLevels,
                   (PosCacheColl**)&(BAMFILE(ext)->pbuffer));
        SEXP res = PROTECT(_pileup_bam(ext, regions, keepFlags,
            reverseComplement, isSimpleCigar, tagFilter, mapqFilter,
            yieldSize, obeyQname,
            asMates, qnamePrefixEnd, qnameSuffixStart, buffer));
        UNPROTECT(2);
        return res;
    }
}
