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

.make_fillers <- function(cigar, ops, letter)
{
    ops_width <- explodeCigarOpLengths(cigar, ops=ops)
    .make_sequence_fillers_from_list_of_widths(ops_width, letter)
}

.D_ranges_on_query_space <- function(cigar, before.hard.clipping=FALSE,
                                            after.soft.clipping=FALSE)
{
    cigarRangesAlongQuerySpace(cigar,
                               before.hard.clipping=before.hard.clipping,
                               after.soft.clipping=after.soft.clipping,
                               ops="D")
}

.N_ranges_on_query_space <- function(cigar, before.hard.clipping=FALSE,
                                            after.soft.clipping=FALSE)
{
    cigarRangesAlongQuerySpace(cigar,
                               before.hard.clipping=before.hard.clipping,
                               after.soft.clipping=after.soft.clipping,
                               ops="N")
}

.D_and_N_ranges_on_query_space <- function(cigar, before.hard.clipping=FALSE,
                                                  after.soft.clipping=FALSE)
{
    D_ranges <- .D_ranges_on_query_space(cigar,
                            before.hard.clipping=before.hard.clipping,
                            after.soft.clipping=after.soft.clipping)
    N_ranges <- .N_ranges_on_query_space(cigar,
                            before.hard.clipping=before.hard.clipping,
                            after.soft.clipping=after.soft.clipping)
    .pcombine(D_ranges, N_ranges)
}

.make_D_and_N_fillers <- function(cigar, D.letter, N.letter)
{
    D_fillers <- .make_fillers(cigar, "D", D.letter)
    N_fillers <- .make_fillers(cigar, "N", N.letter)
    .pcombine(D_fillers, N_fillers)
}

sequenceLayer <- function(x, cigar, from="query", to="reference",
                          D.letter="-", N.letter="-",
                          I.letter="-", S.letter="+", H.letter="+")
{
    if (!is(x, "XStringSet"))
        stop("'x' must be an XStringSet object")
    ## The 8 spaces below are also defined at the top of the src/cigar_utils.c
    ## file in the GenomicRanges package.
    SPACES <- c("reference",
                "reference-N-regions-removed",
                "query",
                "query-before-hard-clipping",
                "query-after-soft-clipping",
                "pairwise",
                "pairwise-N-regions-removed",
                "pairwise-dense")
    from <- match.arg(from, SPACES)
    to <- match.arg(to, SPACES)
    D.letter <- Biostrings:::.normarg_padding.letter(D.letter, seqtype(x))
    N.letter <- Biostrings:::.normarg_padding.letter(N.letter, seqtype(x))
    I.letter <- Biostrings:::.normarg_padding.letter(I.letter, seqtype(x))
    S.letter <- Biostrings:::.normarg_padding.letter(S.letter, seqtype(x))
    H.letter <- Biostrings:::.normarg_padding.letter(H.letter, seqtype(x))
    if (from == to)
        return(x)

    ## Right now, the way 'S.letter' and 'H.letter' are injected in 'x' when
    ## 'to' is "query-before-hard-clipping" can result in padding in the
    ## wrong order (i.e. padding with 'H.letter' followed by padding with
    ## 'S.letter') so we temporarily work around this by enforcing 'S.letter'
    ## and 'H.letter' to be the same.
    if (from != "query" && to == "query-before-hard-clipping" &&
        as.character(S.letter) != as.character(H.letter))
        stop("'H.letter' must be the same as 'S.letter' ",
             "when 'from' is not \"query\" and 'to' ",
             "is \"query-before-hard-clipping\"")

    ## TODO: What follows is a big ugly piece of stinking code (350 lines!).
    ## There is of course a better way...
    inject_at <- ops_to_remove <- NULL
    if (from == "reference") {
        getCigarRanges <- function(cigar, ops)
            cigarRangesAlongReferenceSpace(cigar, ops=ops)
        if (to == "reference-N-regions-removed") {
            ## "reference" -> "reference-N-regions-removed"
            ops_to_remove <- "N"
        } else if (to == "query") {
            ## "reference" -> "query"
            ops_to_remove <- c("D", "N")
            I_inject_at <- getCigarRanges(cigar, "I")
            I_fillers <- .make_fillers(cigar, "I", I.letter)
            S_inject_at <- getCigarRanges(cigar, "S")
            S_fillers <- .make_fillers(cigar, "S", S.letter)
            inject_at <- .pcombine(I_inject_at, S_inject_at)
            fillers <- .pcombine(I_fillers, S_fillers)
        } else if (to == "query-before-hard-clipping") {
            ## "reference" -> "query-before-hard-clipping"
            ops_to_remove <- c("D", "N")
            I_inject_at <- getCigarRanges(cigar, "I")
            I_fillers <- .make_fillers(cigar, "I", I.letter)
            S_inject_at <- getCigarRanges(cigar, "S")
            S_fillers <- .make_fillers(cigar, "S", S.letter)
            H_inject_at <- getCigarRanges(cigar, "H")
            H_fillers <- .make_fillers(cigar, "H", H.letter)
            inject_at <- .pcombine(I_inject_at, S_inject_at)
            inject_at <- .pcombine(inject_at, H_inject_at)
            fillers <- .pcombine(I_fillers, S_fillers)
            fillers <- .pcombine(fillers, H_fillers)
        } else if (to == "query-after-soft-clipping") {
            ## "reference" -> "query-after-soft-clipping"
            ops_to_remove <- c("D", "N")
            ops_to_inject <- "I"
            inject_at <- getCigarRanges(cigar, ops_to_inject)
            fillers <- .make_fillers(cigar, ops_to_inject, I.letter)
        } else if (to == "pairwise") {
            ## "reference" -> "pairwise"
            ops_to_inject <- "I"
            inject_at <- getCigarRanges(cigar, ops_to_inject)
            fillers <- .make_fillers(cigar, ops_to_inject, I.letter)
        } else if (to == "pairwise-N-regions-removed") {
            ## "reference" -> "pairwise-N-regions-removed"
            ops_to_remove <- "N"
            ops_to_inject <- "I"
            inject_at <- getCigarRanges(cigar, ops_to_inject)
            fillers <- .make_fillers(cigar, ops_to_inject, I.letter)
        } else if (to == "pairwise-dense") {
            ## "reference" -> "pairwise-dense"
            ops_to_remove <- c("D", "N")
        }
    } else if (from == "reference-N-regions-removed") {
        getCigarRanges <- function(cigar, ops)
            cigarRangesAlongReferenceSpace(cigar, N.regions.removed=TRUE, ops=ops)
        if (to == "reference") {
            ## "reference-N-regions-removed" -> "reference"
            ops_to_inject <- "N"
            inject_at <- getCigarRanges(cigar, ops_to_inject)
            fillers <- .make_fillers(cigar, ops_to_inject, N.letter)
        } else if (to == "query") {
            ## "reference-N-regions-removed" -> "query"
            ops_to_remove <- "D"
            I_inject_at <- getCigarRanges(cigar, "I")
            I_fillers <- .make_fillers(cigar, "I", I.letter)
            S_inject_at <- getCigarRanges(cigar, "S")
            S_fillers <- .make_fillers(cigar, "S", S.letter)
            inject_at <- .pcombine(I_inject_at, S_inject_at)
            fillers <- .pcombine(I_fillers, S_fillers)
        } else if (to == "query-before-hard-clipping") {
            ## "reference-N-regions-removed" -> "query-before-hard-clipping"
            ops_to_remove <- "D"
            I_inject_at <- getCigarRanges(cigar, "I")
            I_fillers <- .make_fillers(cigar, "I", I.letter)
            S_inject_at <- getCigarRanges(cigar, "S")
            S_fillers <- .make_fillers(cigar, "S", S.letter)
            H_inject_at <- getCigarRanges(cigar, "H")
            H_fillers <- .make_fillers(cigar, "H", H.letter)
            inject_at <- .pcombine(I_inject_at, S_inject_at)
            inject_at <- .pcombine(inject_at, H_inject_at)
            fillers <- .pcombine(I_fillers, S_fillers)
            fillers <- .pcombine(fillers, H_fillers)
        } else if (to == "query-after-soft-clipping") {
            ## "reference-N-regions-removed" -> "query-after-soft-clipping"
            ops_to_remove <- "D"
            ops_to_inject <- "I"
            inject_at <- getCigarRanges(cigar, ops_to_inject)
            fillers <- .make_fillers(cigar, ops_to_inject, I.letter)
        } else if (to == "pairwise") {
            ## "reference-N-regions-removed" -> "pairwise"
            I_inject_at <- getCigarRanges(cigar, "I")
            I_fillers <- .make_fillers(cigar, "I", I.letter)
            N_inject_at <- getCigarRanges(cigar, "N")
            N_fillers <- .make_fillers(cigar, "N", N.letter)
            inject_at <- .pcombine(I_inject_at, N_inject_at)
            fillers <- .pcombine(I_fillers, N_fillers)
        } else if (to == "pairwise-N-regions-removed") {
            ## "reference-N-regions-removed" -> "pairwise-N-regions-removed"
            ops_to_inject <- "I"
            inject_at <- getCigarRanges(cigar, ops_to_inject)
            fillers <- .make_fillers(cigar, ops_to_inject, I.letter)
        } else if (to == "pairwise-dense") {
            ## "reference-N-regions-removed" -> "pairwise-dense"
            ops_to_remove <- "D"
        }
    } else if (from == "query") {
        getCigarRanges <- function(cigar, ops)
            cigarRangesAlongQuerySpace(cigar, ops=ops)
        if (to == "reference") {
            ## "query" -> "reference"
            ops_to_remove <- c("I", "S")
            inject_at <- .D_and_N_ranges_on_query_space(cigar)
            fillers <- .make_D_and_N_fillers(cigar, D.letter, N.letter)
        } else if (to == "reference-N-regions-removed") {
            ## "query" -> "reference-N-regions-removed"
            ops_to_remove <- c("S", "I")
            ops_to_inject <- "D"
            inject_at <- getCigarRanges(cigar, ops_to_inject)
            fillers <- .make_fillers(cigar, ops_to_inject, D.letter)
        } else if (to == "query-before-hard-clipping") {
            ## "query" -> "query-before-hard-clipping"
            ops_to_inject <- "H"
            inject_at <- getCigarRanges(cigar, ops_to_inject)
            fillers <- .make_fillers(cigar, ops_to_inject, H.letter)
        } else if (to == "query-after-soft-clipping") {
            ## "query" -> "query-after-soft-clipping"
            ops_to_remove <- "S"
        } else if (to == "pairwise") {
            ## "query" -> "pairwise"
            ops_to_remove <- "S"
            inject_at <- .D_and_N_ranges_on_query_space(cigar)
            fillers <- .make_D_and_N_fillers(cigar, D.letter, N.letter)
        } else if (to == "pairwise-N-regions-removed") {
            ## "query" -> "pairwise-N-regions-removed"
            ops_to_remove <- "S"
            ops_to_inject <- "D"
            inject_at <- getCigarRanges(cigar, ops_to_inject)
            fillers <- .make_fillers(cigar, ops_to_inject, D.letter)
        } else if (to == "pairwise-dense") {
            ## "query" -> "pairwise-dense"
            ops_to_remove <- c("I", "S")
        }
    } else if (from == "query-before-hard-clipping") {
        getCigarRanges <- function(cigar, ops)
            cigarRangesAlongQuerySpace(cigar, before.hard.clipping=TRUE, ops=ops)
        if (to == "reference") {
            ## "query-before-hard-clipping" -> "reference"
            ops_to_remove <- c("H", "S", "I")
            inject_at <- .D_and_N_ranges_on_query_space(cigar,
                                                 before.hard.clipping=TRUE)
            fillers <- .make_D_and_N_fillers(cigar, D.letter, N.letter)
        } else if (to == "reference-N-regions-removed") {
            ## "query-before-hard-clipping" -> "reference-N-regions-removed"
            ops_to_remove <- c("H", "S", "I")
            ops_to_inject <- "D"
            inject_at <- getCigarRanges(cigar, ops_to_inject)
            fillers <- .make_fillers(cigar, ops_to_inject, D.letter)
        } else if (to == "query") {
            ## "query-before-hard-clipping" -> "query"
            ops_to_remove <- "H"
        } else if (to == "query-after-soft-clipping") {
            ## "query-before-hard-clipping" -> "query-after-soft-clipping"
            ops_to_remove <- c("H", "S")
        } else if (to == "pairwise") {
            ## "query-before-hard-clipping" -> "pairwise"
            ops_to_remove <- c("H", "S")
            inject_at <- .D_and_N_ranges_on_query_space(cigar,
                                                 before.hard.clipping=TRUE)
            fillers <- .make_D_and_N_fillers(cigar, D.letter, N.letter)
        } else if (to == "pairwise-N-regions-removed") {
            ## "query-before-hard-clipping" -> "pairwise-N-regions-removed"
            ops_to_remove <- c("H", "S")
            ops_to_inject <- "D"
            inject_at <- getCigarRanges(cigar, ops_to_inject)
            fillers <- .make_fillers(cigar, ops_to_inject, D.letter)
        } else if (to == "pairwise-dense") {
            ## "query-before-hard-clipping" -> "pairwise-dense"
            ops_to_remove <- c("I", "S", "H")
        }
    } else if (from == "query-after-soft-clipping") {
        getCigarRanges <- function(cigar, ops)
            cigarRangesAlongQuerySpace(cigar, after.soft.clipping=TRUE, ops=ops)
        if (to == "reference") {
            ## "query-after-soft-clipping" -> "reference"
            ops_to_remove <- "I"
            inject_at <- .D_and_N_ranges_on_query_space(cigar,
                                                 after.soft.clipping=TRUE)
            fillers <- .make_D_and_N_fillers(cigar, D.letter, N.letter)
        } else if (to == "reference-N-regions-removed") {
            ## "query-after-soft-clipping" -> "reference-N-regions-removed"
            ops_to_remove <- "I"
            ops_to_inject <- "D"
            inject_at <- getCigarRanges(cigar, ops_to_inject)
            fillers <- .make_fillers(cigar, ops_to_inject, D.letter)
        } else if (to == "query") {
            ## "query-after-soft-clipping" -> "query"
            ops_to_inject <- "S"
            inject_at <- getCigarRanges(cigar, ops_to_inject)
            fillers <- .make_fillers(cigar, ops_to_inject, S.letter)
        } else if (to == "query-before-hard-clipping") {
            ## "query-after-soft-clipping" -> "query-before-hard-clipping"
            S_inject_at <- getCigarRanges(cigar, "S")
            S_fillers <- .make_fillers(cigar, "S", S.letter)
            H_inject_at <- getCigarRanges(cigar, "H")
            H_fillers <- .make_fillers(cigar, "H", H.letter)
            inject_at <- .pcombine(S_inject_at, H_inject_at)
            fillers <- .pcombine(S_fillers, H_fillers)
        } else if (to == "pairwise") {
            ## "query-after-soft-clipping" -> "pairwise"
            inject_at <- .D_and_N_ranges_on_query_space(cigar,
                                                 after.soft.clipping=TRUE)
            fillers <- .make_D_and_N_fillers(cigar, D.letter, N.letter)
        } else if (to == "pairwise-N-regions-removed") {
            ## "query-after-soft-clipping" -> "pairwise-N-regions-removed"
            ops_to_inject <- "D"
            inject_at <- getCigarRanges(cigar, ops_to_inject)
            fillers <- .make_fillers(cigar, ops_to_inject, D.letter)
        } else if (to == "pairwise-dense") {
            ## "query-after-soft-clipping" -> "pairwise-dense"
            ops_to_remove <- "I"
        }
    } else if (from == "pairwise") {
        getCigarRanges <- function(cigar, ops)
            cigarRangesAlongPairwiseSpace(cigar, ops=ops)
        if (to == "reference") {
            ## "pairwise" -> "reference"
            ops_to_remove <- "I"
        } else if (to == "reference-N-regions-removed") {
            ## "pairwise" -> "reference-N-regions-removed"
            ops_to_remove <- c("I", "N")
        } else if (to == "query") {
            ## "pairwise" -> "query"
            ops_to_remove <- c("D", "N")
            ops_to_inject <- "S"
            inject_at <- getCigarRanges(cigar, ops_to_inject)
            fillers <- .make_fillers(cigar, ops_to_inject, S.letter)
        } else if (to == "query-before-hard-clipping") {
            ## "pairwise" -> "query-before-hard-clipping"
            ops_to_remove <- c("D", "N")
            S_inject_at <- getCigarRanges(cigar, "S")
            S_fillers <- .make_fillers(cigar, "S", S.letter)
            H_inject_at <- getCigarRanges(cigar, "H")
            H_fillers <- .make_fillers(cigar, "H", H.letter)
            inject_at <- .pcombine(S_inject_at, H_inject_at)
            fillers <- .pcombine(S_fillers, H_fillers)
        } else if (to == "query-after-soft-clipping") {
            ## "pairwise" -> "query-after-soft-clipping"
            ops_to_remove <- c("D", "N")
        } else if (to == "pairwise-N-regions-removed") {
            ## "pairwise" -> "pairwise-N-regions-removed"
            ops_to_remove <- "N"
        } else if (to == "pairwise-dense") {
            ## "pairwise" -> "pairwise-dense"
            ops_to_remove <- c("I", "D", "N")
        }
    } else if (from == "pairwise-N-regions-removed") {
        getCigarRanges <- function(cigar, ops)
            cigarRangesAlongPairwiseSpace(cigar, N.regions.removed=TRUE, ops=ops)
        if (to == "reference") {
            ## "pairwise-N-regions-removed" -> "reference"
            ops_to_remove <- "I"
            ops_to_inject <- "N"
            inject_at <- getCigarRanges(cigar, ops_to_inject)
            fillers <- .make_fillers(cigar, ops_to_inject, N.letter)
        } else if (to == "reference-N-regions-removed") {
            ## "pairwise-N-regions-removed" -> "reference-N-regions-removed"
            ops_to_remove <- "I"
        } else if (to == "query") {
            ## "pairwise-N-regions-removed" -> "query"
            ops_to_remove <- "D"
            ops_to_inject <- "S"
            inject_at <- getCigarRanges(cigar, ops_to_inject)
            fillers <- .make_fillers(cigar, ops_to_inject, S.letter)
        } else if (to == "query-before-hard-clipping") {
            ## "pairwise-N-regions-removed" -> "query-before-hard-clipping"
            ops_to_remove <- "D"
            S_inject_at <- getCigarRanges(cigar, "S")
            S_fillers <- .make_fillers(cigar, "S", S.letter)
            H_inject_at <- getCigarRanges(cigar, "H")
            H_fillers <- .make_fillers(cigar, "H", H.letter)
            inject_at <- .pcombine(S_inject_at, H_inject_at)
            fillers <- .pcombine(S_fillers, H_fillers)
        } else if (to == "query-after-soft-clipping") {
            ## "pairwise-N-regions-removed" -> "query-after-soft-clipping"
            ops_to_remove <- "D"
        } else if (to == "pairwise") {
            ## "pairwise-N-regions-removed" -> "pairwise"
            ops_to_inject <- "N"
            inject_at <- getCigarRanges(cigar, ops_to_inject)
            fillers <- .make_fillers(cigar, ops_to_inject, N.letter)
        } else if (to == "pairwise-dense") {
            ## "pairwise-N-regions-removed" -> "pairwise-dense"
            ops_to_remove <- c("I", "D")
        }
    } else if (from == "pairwise-dense") {
        getCigarRanges <- function(cigar, ops)
            cigarRangesAlongPairwiseSpace(cigar, dense=TRUE, ops=ops)
        if (to == "reference") {
            ## "pairwise-dense" -> "reference"
            D_inject_at <- getCigarRanges(cigar, "D")
            D_fillers <- .make_fillers(cigar, "D", D.letter)
            N_inject_at <- getCigarRanges(cigar, "N")
            N_fillers <- .make_fillers(cigar, "N", N.letter)
            inject_at <- .pcombine(D_inject_at, N_inject_at)
            fillers <- .pcombine(D_fillers, N_fillers)
        } else if (to == "reference-N-regions-removed") {
            ## "pairwise-dense" -> "reference-N-regions-removed"
            ops_to_inject <- "D"
            inject_at <- getCigarRanges(cigar, ops_to_inject)
            fillers <- .make_fillers(cigar, ops_to_inject, D.letter)
        } else if (to == "query") {
            ## "pairwise-dense" -> "query"
            I_inject_at <- getCigarRanges(cigar, "I")
            I_fillers <- .make_fillers(cigar, "I", I.letter)
            S_inject_at <- getCigarRanges(cigar, "S")
            S_fillers <- .make_fillers(cigar, "S", S.letter)
            inject_at <- .pcombine(I_inject_at, S_inject_at)
            fillers <- .pcombine(I_fillers, S_fillers)
        } else if (to == "query-before-hard-clipping") {
            ## "pairwise-dense" -> "query-before-hard-clipping"
            I_inject_at <- getCigarRanges(cigar, "I")
            I_fillers <- .make_fillers(cigar, "I", I.letter)
            S_inject_at <- getCigarRanges(cigar, "S")
            S_fillers <- .make_fillers(cigar, "S", S.letter)
            H_inject_at <- getCigarRanges(cigar, "H")
            H_fillers <- .make_fillers(cigar, "H", H.letter)
            inject_at <- .pcombine(.pcombine(I_inject_at, S_inject_at),
                                   H_inject_at)
            fillers <- .pcombine(.pcombine(I_fillers, S_fillers), H_fillers)
        } else if (to == "query-after-soft-clipping") {
            ## "pairwise-dense" -> "query-after-soft-clipping"
            ops_to_inject <- "I"
            inject_at <- getCigarRanges(cigar, ops_to_inject)
            fillers <- .make_fillers(cigar, ops_to_inject, I.letter)
        } else if (to == "pairwise") {
            ## "pairwise-dense" -> "pairwise"
            I_inject_at <- getCigarRanges(cigar, "I")
            I_fillers <- .make_fillers(cigar, "I", I.letter)
            D_inject_at <- getCigarRanges(cigar, "D")
            D_fillers <- .make_fillers(cigar, "D", D.letter)
            N_inject_at <- getCigarRanges(cigar, "N")
            N_fillers <- .make_fillers(cigar, "N", N.letter)
            inject_at <- .pcombine(.pcombine(I_inject_at, D_inject_at),
                                   N_inject_at)
            fillers <- .pcombine(.pcombine(I_fillers, D_fillers), N_fillers)
        } else if (to == "pairwise-N-regions-removed") {
            ## "pairwise-dense" -> "pairwise-N-regions-removed"
            I_inject_at <- getCigarRanges(cigar, "I")
            I_fillers <- .make_fillers(cigar, "I", I.letter)
            D_inject_at <- getCigarRanges(cigar, "D")
            D_fillers <- .make_fillers(cigar, "D", D.letter)
            inject_at <- .pcombine(I_inject_at, D_inject_at)
            fillers <- .pcombine(I_fillers, D_fillers)
        }
    }

    at <- inject_at
    if (is.null(inject_at)) {
        value <- NULL
    } else {
        value <- fillers
    }
    if (length(ops_to_remove) != 0L) {
        at2 <- getCigarRanges(cigar, ops_to_remove)
        value2 <- .make_empty_sequences(at2, class=class(x))
        at <- .pcombine(at, at2)
        value <- .pcombine(value, value2)
    }
    replaceAt(x, at, value=value)
}

