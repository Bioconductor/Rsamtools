### =========================================================================
### stackStringsFromBam() and related
### -------------------------------------------------------------------------


.make_empty_sequences <- function(skeleton, class="BStringSet")
{
    if (is.null(skeleton))
        return(NULL)
    skeleton <- PartitioningByEnd(skeleton)
    skeleton_len <- length(skeleton)
    if (skeleton_len == 0L) {
        unlisted_len <- 0L
    } else {
        unlisted_len <- end(skeleton)[skeleton_len]
    }
    unlisted_ans <- rep.int(as("", class), unlisted_len)
    relist(unlisted_ans, skeleton)
}

### 'filler_width' must be an integer vector, and 'letter' an XString object
### of length 1.
.make_sequence_fillers_from_widths <- function(filler_width, letter)
{
    if (length(filler_width) == 0L) {
        max_width <- 0L
        at <- IRanges()
    } else {
        max_width <- max(filler_width)
        at <- IRanges(1L, filler_width)
    }
    biggest_filler <- rep.int(letter, max_width)
    extractAt(biggest_filler, at)
}

### 'filler_widths' must be an IntegerList object (or list of integers).
.make_sequence_fillers_from_list_of_widths <- function(filler_widths, letter)
{
    unlisted_widths <- unlist(filler_widths, use.names=FALSE)
    unlisted_ans <- .make_sequence_fillers_from_widths(unlisted_widths, letter)
    relist(unlisted_ans, filler_widths)
}

### Parallel combine.
.pcombine <- function(x, y)
{
    if (is.null(x))
        return(y)
    if (is.null(y))
        return(x)
    if (length(x) != length(y))
        stop("'x' and 'y' must have the same length")
    xy <- c(x, y)
    collate_subscript <- IRanges:::make_XYZxyz_to_XxYyZz_subscript(length(x))
    ans_flesh <- unlist(xy[collate_subscript], use.names=FALSE)
    ans_breakpoints <- end(PartitioningByEnd(x)) + end(PartitioningByEnd(y))
    ans_skeleton <- PartitioningByEnd(ans_breakpoints)
    relist(ans_flesh, ans_skeleton)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### sequenceLayer()
###

.D_and_N_ranges_on_query_space <- function(cigar)
{
    ans1 <- cigarRangesOnQuerySpace(cigar, ops="D")
    ans2 <- cigarRangesOnQuerySpace(cigar, ops="N")
    .pcombine(ans1, ans2)
}

.D_and_N_fillers <- function(cigar, D.letter, N.letter)
{
    D_width <- width(cigarRangesOnReferenceSpace(cigar, ops="D"))
    N_width <- width(cigarRangesOnReferenceSpace(cigar, ops="N"))
    ans1 <- .make_sequence_fillers_from_list_of_widths(D_width, D.letter)
    ans2 <- .make_sequence_fillers_from_list_of_widths(N_width, N.letter)
    .pcombine(ans1, ans2)
}

.I_or_S_fillers <- function(cigar, ops, letter)
{
    ops_width <- width(cigarRangesOnQuerySpace(cigar, ops=ops))
    .make_sequence_fillers_from_list_of_widths(ops_width, letter)
}

sequenceLayer <- function(x, cigar, from="query", to="reference",
                          D.letter="-", N.letter="-")
{
    if (!is(x, "XStringSet"))
        stop("'x' must be an XStringSet object")
    LAYOUTS <- c("reference", "query", "pairwise")
    from <- match.arg(from, LAYOUTS)
    to <- match.arg(to, LAYOUTS)
    D.letter <- Biostrings:::.normarg_padding.letter(D.letter, seqtype(x))
    N.letter <- Biostrings:::.normarg_padding.letter(N.letter, seqtype(x))
    if (from == to)
        return(x)
    if (from == "query") {
        if (to == "reference") {
            ## From "query" to "reference".
            ops_to_remove <- c("I", "S")
        } else {
            ## From "query" to "pairwise".
            ops_to_remove <- "S"
        }
        remove_at <- cigarRangesOnQuerySpace(cigar, ops=ops_to_remove)
        inject_at <- .D_and_N_ranges_on_query_space(cigar)
        fillers <- .D_and_N_fillers(cigar, D.letter, N.letter)
    } else if (from == "reference") {
        if (to == "query") {
            ## From "reference" to "query".
            ops_to_remove <- c("D", "N")
            remove_at <- cigarRangesOnReferenceSpace(cigar, ops=ops_to_remove)
            ops_to_inject <- c("I", "S")
            inject_at <- cigarRangesOnReferenceSpace(cigar, ops=ops_to_inject)
            letter <- as("-", class(x))[[1L]]
            fillers <- .I_or_S_fillers(cigar, ops_to_inject, letter)
        } else {
            ## From "reference" to "pairwise".
            remove_at <- NULL
            ops_to_inject <- "I"
            inject_at <- cigarRangesOnReferenceSpace(cigar, ops=ops_to_inject)
            letter <- as("-", class(x))[[1L]]
            fillers <- .I_or_S_fillers(cigar, ops_to_inject, letter)
        }
    } else {
        if (to == "query") {
            ## From "pairwise" to "query".
            ops_to_remove <- c("D", "N")
            remove_at <- cigarRangesOnPairwiseSpace(cigar, ops=ops_to_remove)
            ops_to_inject <- "S"
            inject_at <- cigarRangesOnPairwiseSpace(cigar, ops=ops_to_inject)
            letter <- as("-", class(x))[[1L]]
            fillers <- .I_or_S_fillers(cigar, ops_to_inject, letter)
        } else {
            ## From "pairwise" to "reference".
            ops_to_remove <- "I"
            remove_at <- cigarRangesOnPairwiseSpace(cigar, ops=ops_to_remove)
            inject_at <- fillers <- NULL
        }
    }
    at <- .pcombine(remove_at, inject_at)
    empty_sequences <- .make_empty_sequences(remove_at, class=class(x))
    value <- .pcombine(empty_sequences, fillers)
    replaceAt(x, at, value=value)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### stackStringsFromBam()
###

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

if (FALSE) {

library(Rsamtools)

qseq <- DNAStringSet(c(seq1="nnAAAACCTTGAAA", seq2="AAAAATCCCTTTTnn"))
cigar <- c("2S4M1D2M2I1M5N3M", "5M1I3M2D4M2S")

### Check compatibility between 'qseq' and 'cigar' (because this is not cheap,
### sequenceLayer() doesn't do it):
stopifnot(identical(width(qseq), cigarToQWidth(cigar)))

rseq0 <- DNAStringSet(c(seq1="AAAA-CCG-----AAA", seq2="AAAAACCC--TTTT"))
pseq0 <- DNAStringSet(c(seq1="AAAA-CCTTG-----AAA", seq2="AAAAATCCC--TTTT"))
### 1 quick sanity check before we start testing sequenceLayer():
stopifnot(identical(width(rseq0), cigarToWidth(cigar)))

###              ~ ~ ~ ~ ~ Testing sequenceLayer() ~ ~ ~ ~ ~
### query <-> reference
rseq <- sequenceLayer(qseq, cigar, from="query", to="reference")
stopifnot(all(rseq0 == rseq))
qseq2 <- sequenceLayer(rseq, cigar, from="reference", to="query")
stopifnot(all(width(qseq) == width(qseq2)))
rseq2 <- sequenceLayer(qseq2, cigar, from="query", to="reference")
stopifnot(all(rseq == rseq2))
### query <-> pairwise
pseq <- sequenceLayer(qseq, cigar, from="query", to="pairwise")
stopifnot(all(pseq0 == pseq))
qseq3 <- sequenceLayer(pseq, cigar, from="pairwise", to="query")
stopifnot(all(width(qseq) == width(qseq3)))
pseq2 <- sequenceLayer(qseq3, cigar, from="query", to="pairwise")
stopifnot(all(pseq == pseq2))
### lay query and reference to pairwise space and compare
pseq3a <- sequenceLayer(qseq2, cigar, from="query", to="pairwise")
pseq3b <- sequenceLayer(rseq, cigar, from="reference", to="pairwise")
stopifnot(all(pseq3a == pseq3b))


bamfile <- BamFile(system.file("extdata", "ex1.bam", package="Rsamtools"))

stackStringsFromBam(bamfile, GRanges("seq1", IRanges(1, 60)))

options(showHeadLines=25)
options(showTailLines=2)
stackStringsFromBam(bamfile, GRanges("seq1", IRanges(61, 120)))

stacked_reads <- stackStringsFromBam(bamfile, "seq2:1509-1519")
stacked_reads  # deletion in read 13
stackStringsFromBam(bamfile, "seq2:1509-1519", what="qual")
consensusMatrix(stacked_reads)


library(RNAseqData.HNRNPC.bam.chr14)
bamfile <- BamFile(RNAseqData.HNRNPC.bam.chr14_BAMFILES[1])

my_ROI <- GRanges("chr14", IRanges(19650095, 19650159)) # my Region Of Interest
readGAlignments(bamfile, param=ScanBamParam(which=my_ROI))
stackStringsFromBam(bamfile, my_ROI)
}

