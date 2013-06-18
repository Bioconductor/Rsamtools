### =========================================================================
### sequenceLayer()
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
}

