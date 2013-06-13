### =========================================================================
### stackStringsFromBam() and related
### -------------------------------------------------------------------------


.make_empty_sequences <- function(skeleton, class="BStringSet")
{
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

### 'filler_width' must be an integer vector.
.make_sequence_fillers_from_widths <- function(filler_width, letter,
                                               class="BStringSet")
{
    if (length(filler_width) == 0L) {
        max_width <- 0L
        at <- IRanges()
    } else {
        max_width <- max(filler_width)
        at <- IRanges(1L, filler_width)
    }
    biggest_filler <- rep.int(as(letter, class)[[1L]], max_width)
    extractAt(biggest_filler, at)
}

### 'filler_widths' must be an IntegerList object (or list of integers).
.make_sequence_fillers_from_list_of_widths <- function(filler_widths, letter,
                                                       class="BStringSet")
{
    unlisted_widths <- unlist(filler_widths, use.names=FALSE)
    unlisted_ans <- .make_sequence_fillers_from_widths(unlisted_widths, letter,
                                                       class=class)
    relist(unlisted_ans, filler_widths)
}

### Parallel combine.
.pcombine <- function(x, y)
{
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
### TODO: (1) Add option for dropping sequences that are not visible.
###       (2) Let the user choose the letter to use for filling deletions
###           and gaps.
###       (3) If the letter for padding and gaps is the same (e.g. "+"), then
###           doing consensusMatrix() on the result of sequenceLayer() should
###           give a result consistent with coverage().
###

sequenceLayer <- function(x, cigar, layout="query-to-reference")
{
    if (!is(x, "XStringSet"))
        stop("'x' must be an XStringSet object")
    LAYOUTS <- c("query-to-reference",
                 "reference-to-query",
                 "aligned-query",
                 "aligned-reference")
    layout <- match.arg(layout, LAYOUTS)

    ops1 <- c("I", "S")
    ops2 <- c("D", "N")
    if (layout == "query-to-reference") {
        del_ranges <- cigarRangesOnQuery(cigar, ops=ops1)
        fill_ranges <- cigarRangesOnQuery(cigar, ops=ops2)
        to_ranges <- cigarRangesOnReference(cigar,
                                            pos=rep.int(1L, length(cigar)),
                                            ops=ops2)
    } else if (layout == "reference-to-query") {
        del_ranges <- cigarRangesOnReference(cigar,
                                             pos=rep.int(1L, length(cigar)),
                                             ops=ops2)
        fill_ranges <- cigarRangesOnReference(cigar,
                                              pos=rep.int(1L, length(cigar)),
                                              ops=ops1)
        to_ranges <- cigarRangesOnQuery(cigar, ops=ops1)
    } else {
        stop("\"", layout, "\" layout is not implemented yet, sorry!")
    }
    at <- .pcombine(del_ranges, fill_ranges)
    empty_sequences <- .make_empty_sequences(del_ranges, class=class(x))
    fillers <- .make_sequence_fillers_from_list_of_widths(width(to_ranges), "-",
                                                          class=class(x))
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

stackStringsFromBam <- function(file, param, use.names=FALSE,
                                what="seq", padding.letter=NA)
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
    if (what == "seq") {
        if (identical(padding.letter, NA))
            padding.letter <- DNAString("+")
    } else {  # what == "qual"
        what_col <- BStringSet(what_col)
        if (identical(padding.letter, NA))
            padding.letter <- BString(" ")
    }
    layed_seq <- sequenceLayer(what_col, cigar(gal))
    ans <- stackStrings(layed_seq, start(region), end(region), padding.letter,
                        shift=start(gal)-1L)
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

x <- DNAStringSet(c(seq1="AAAA", seq2="AAAAATCCCTTTTNN"))
cigar <- c("4M", "5M1I3M2D4M2S")

### Check compatibility between 'x' and 'cigar' (because this is not cheap,
### sequenceLayer() doesn't do it):
stopifnot(identical(width(x), cigarToQWidth(cigar)))

target <- DNAStringSet(c(seq1="AAAA", seq2="AAAAACCC--TTTT"))
### A quick sanity check before we start testing sequenceLayer():
stopifnot(identical(width(target), cigarToWidth(cigar)))

### Testing sequenceLayer():
x2 <- sequenceLayer(x, cigar, layout="query-to-reference")
stopifnot(all(target == x2))
x3 <- sequenceLayer(x2, cigar, layout="reference-to-query")
stopifnot(all(width(x) == width(x3)))
x4 <- sequenceLayer(x3, cigar, layout="query-to-reference")
stopifnot(all(x2 == x4))



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

