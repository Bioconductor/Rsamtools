### =========================================================================
### stackSequences() and related
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
###       (2) Let the user choose the letter used for padding and for filling
###           deletions and gaps.
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
### stackSequences()
###

setGeneric("stackSequences", signature="x",
    function(x, from, to, pos=1L, pad.letter="+")
        standardGeneric("stackSequences")
)

setMethod("stackSequences", "XStringSet",
    function(x, from, to, pos=1L, pad.letter="+")
    {
        if (!isSingleNumber(from))
            stop("'from' must be a single integer")
        if (!is.integer(from))
            from <- as.integer(from)
        if (!isSingleNumber(to))
            stop("'to' must be a single integer")
        if (!is.integer(to))
            to <- as.integer(to)
        width0 <- to - from + 1L
        if (width0 < 0L)
            stop("'to' must be >= 'from - 1L'")
        if (!is.numeric(pos))
            stop("'pos' must be a vector of integers")
        if (!is.integer(pos)) 
            pos <- as.integer(pos)
        pos <- Biostrings:::.V_recycle(pos, x, "pos", "'length(x)'")
        if (!isSingleString(pad.letter) || nchar(pad.letter) != 1L)
            stop("'pad.letter' must be a single letter")

        left_margin <- pos - from
        right_margin <- to - (pos + width(x) - 1L)

        left_pad <- pmin(pmax(left_margin, 0L), width0)
        right_pad <- pmin(pmax(right_margin, 0L), width0)
        left_trim <- pmin(pmax(-left_margin, 0L), width(x))
        right_trim <- pmin(pmax(-right_margin, 0L), width(x))

        left <- .make_sequence_fillers_from_widths(left_pad, pad.letter,
                                                   class=class(x))
        right <- .make_sequence_fillers_from_widths(right_pad, pad.letter,
                                                    class=class(x))
        middle <- narrow(x, start=1L+left_trim, end=-(1L+right_trim))
        xscat(left, middle, right)
    }
)

setMethod("stackSequences", "BamFile",
    function(x, from, to, pos=1L, pad.letter="+")
    {
        if (!is(from, "GenomicRanges"))
            stop("'from' must be a GRanges object ",
                 "when 'x' is a BamFile object")
        if (!missing(to))
            warning("'to' is ignored when 'x' is a BamFile object")
        if (!identical(pos, 1L))
            warning("'pos' is ignored when 'x' is a BamFile object")
        param <- ScanBamParam(what="seq", which=from)
        gal <- readGAlignmentsFromBam(x, param=param)
        layed_seq <- sequenceLayer(mcols(gal)$seq, cigar(gal))
        stackSequences(layed_seq,
                       from=start(from), to=end(from), pos=start(gal))
    }
)

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

stackSequences(bamfile, GRanges("seq1", IRanges(1, 60)))

options(showHeadLines=25)
options(showTailLines=2)
stackSequences(bamfile, GRanges("seq1", IRanges(61, 120)))

stacked_reads <- stackSequences(bamfile, GRanges("seq2", IRanges(1509, 1519)))
stacked_reads  # deletion in read 13
consensusMatrix(stacked_reads)


library(RNAseqData.HNRNPC.bam.chr14)
bamfile <- BamFile(RNAseqData.HNRNPC.bam.chr14_BAMFILES[1])

my_ROI <- GRanges("chr14", IRanges(19650095, 19650159)) # my Region Of Interest
readGAlignments(bamfile, param=ScanBamParam(which=my_ROI))
stackSequences(bamfile, my_ROI)
}

