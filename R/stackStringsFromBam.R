### =========================================================================
### stackStringsFromBam()
### -------------------------------------------------------------------------


### Should always return a ScanBamParam object containing exactly 1 genomic
### region.
.normarg_param <- function(param)
{
    if (isSingleString(param)) {
        tmp1 <- strsplit(param, ":", fixed=TRUE)[[1L]]
        if (length(tmp1) != 2L)
            stop("when a character string, 'param' must be ",
                 "of the form \"chr14:5201-5300\"")
        tmp2 <- as.integer(strsplit(tmp1[2L], "-", fixed=TRUE)[[1L]])
        if (length(tmp2) != 2L || any(is.na(tmp2)))
            stop("when a character string, 'param' must be ",
                 "of the form \"chr14:5201-5300\"")
        param <- GRanges(tmp1[1L], IRanges(tmp2[1L], tmp2[2L]))
    }
    if (is(param, "GenomicRanges")) {
        if (length(param) != 1L)
            stop("when a GRanges object, 'param' must have length 1")
        param <- ScanBamParam(which=param)
        return(param)
    }
    if (is(param, "RangesList")) {
        ## We support RangesList just because ScanBamParam() supports it too
        ## and also because that's what's returned by bamWhich().
        region <- unlist(param, use.names=FALSE)
        if (length(region) != 1L)
            stop("when a RangesList object, 'param' must contain exactly 1 ",
                 "genomic region\n  (i.e. 'unlist(param)' must have length 1)")
        param <- ScanBamParam(which=param)
        return(param)
    }
    if (!is(param, "ScanBamParam"))
        stop("'param' must be either a ScanBamParam or RangesList object ",
             "containing\n  exactly 1 genomic region, or a GRanges object ",
             "of length 1, or a character\n  string specifying a single ",
             "genomic region (in the \"chr14:5201-5300\" format)")
    region <- unlist(bamWhich(param), use.names=FALSE)
    if (length(region) != 1L)
        stop("when a ScanBamParam object, 'param' must contain exactly 1 ",
             "genomic region\n  (i.e. 'unlist(bamWhich(param))' must have ",
             "length 1)")
    param
}

### NOTE: Using the same letter ("-") for gaps (N) and deletions (D) is
###       problematic when one tries to infer coverage by calling
###       consensusMatrix() on the result of sequenceLayer(): then gaps
###       and deletions will both generate coverage when only the latters
###       should! Hence the result won't be consistent with coverage().
###       OTOH using "+" for gaps is compatible with coverage() but is
###       not visually appealing (the gaps, which are usually big, cannot
###       be seen anymore).
###       "-" and "+" are the only 2 special letters in Biostrings DNA
###       alphabet (DNA_ALPHABET). Should we add 1 more? Should it be " ",
###       ".", or "*"?
stackStringsFromBam <- function(file, param, use.names=FALSE, what="seq",
                                D.letter="-", N.letter="-",
                                Lpadding.letter="+", Rpadding.letter="+")
{
    param <- .normarg_param(param)
    region <- unlist(bamWhich(param), use.names=FALSE)
    what <- match.arg(what, c("seq", "qual"))
    param_what <- bamWhat(param)
    if (!(what %in% param_what))
        bamWhat(param) <- c(param_what, what)
    gal <- readGAlignmentsFromBam(file, use.names=use.names, param=param)
    gal_mcols <- mcols(gal)
    what_col_idx <- match(what, colnames(gal_mcols))
    what_col <- gal_mcols[[what_col_idx]]
    if (what == "qual")
        what_col <- BStringSet(what_col)
    layed_seq <- sequenceLayer(what_col, cigar(gal),
                               D.letter=D.letter, N.letter=N.letter)
    ans <- stackStrings(layed_seq, start(region), end(region),
                        shift=start(gal)-1L,
                        Lpadding.letter=Lpadding.letter,
                        Rpadding.letter=Rpadding.letter)
    if (!(what %in% param_what)) {
        ## Remove the what column from 'gal_mcols'.
        gal_mcols <- gal_mcols[ , -what_col_idx, drop=FALSE]
        ## Sadly, subsetting a DataFrame will mangle the colnames of the
        ## returned DataFrame if it has duplicated colnames. Since we of
        ## course don't want this, we fix them.
        colnames(gal_mcols) <- param_what
    }
    names(ans) <- names(gal)
    mcols(ans) <- gal_mcols
    ans
}

