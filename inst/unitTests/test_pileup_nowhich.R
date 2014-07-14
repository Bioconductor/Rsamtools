suppressMessages({
    library(Rsamtools)
    library(RUnit)
})

.reorder_data.frame <- function(df) {
    df <- df[do.call(order,df),]
    rownames(df) <- NULL
    df
}

.unordered_check <- function(expected, xx) {
    xx <- .reorder_data.frame(xx)
    expected <- .reorder_data.frame(expected)
    checkIdentical(xx, expected)
}

.s_levels <- function() {
    levels(strand())
}
.n_levels <- function() {
    c("A", "C", "G", "T", "N", "=", "-")
}

.tiny.sam_seqlevels <- function() {
    fl <- system.file(package="Rsamtools", "extdata", "tiny.bam")
    bf <- BamFile(fl)
    seqlevels(bf)
}
.ex1.sam_seqlevels <- function() {
    fl <- system.file(package="Rsamtools", "extdata", "ex1.bam")
    bf <- BamFile(fl)
    seqlevels(bf)
}
.no_which.sam_seqlevels <- function() {
    fl <- system.file(package="Rsamtools", "extdata", "no_which_whole_file.bam")
    bf <- BamFile(fl)
    seqlevels(bf)
}

## target as data.frame
.tadf <- function(pos=NULL, strand=NULL, nucleotide=NULL, cycle_bin=NULL,
                  count=NULL, which_label=NULL, seqnames=NULL,
                  use_ex1.bam_levels=FALSE) {
    if(is.null(pos) || is.null(count) || is.null(seqnames))
        stop("'pos', 'count', and 'seqnames' must not be 'NULL'")
    target <- data.frame(pos=as.integer(pos), stringsAsFactors=FALSE)
    seqnames_levels <- .no_which.sam_seqlevels()
    target <- cbind(seqnames=factor(seqnames, levels=seqnames_levels), target)
    
    if(!is.null(strand))
        target <- cbind(target, strand=factor(strand, levels=.s_levels()))
    if(!is.null(nucleotide))
        target <- cbind(target,
                        nucleotide=factor(nucleotide, levels=.n_levels()))
    if(!is.null(cycle_bin))
        target <- cbind(target,
                        cycle_bin=cycle_bin)
    target <- cbind(target, count=as.integer(count))
    if(!is.null(which_label))
        target <- cbind(target, which_label=which_label)
    target
}

## make which labels
.mwls <- function(param, run_lens) {
    wl <- Rsamtools:::.scanBam_extract_which_labels(param)
    if(length(wl) != length(run_lens))
        stop("mismatched lengths, length(wl): '%d' length(run_lens): '%d'\n",
             length(wl), length(run_lens))
    if(all(run_lens == 0L))
        factor(levels=wl)
    else
        rep(wl, run_lens)
}

## END HELPER FUNCTIONS

## file helpers
no_which_whole_file <-
    system.file(package="Rsamtools", "extdata", "no_which_whole_file.bam")

no_which_buf_pileup <-
    system.file(package="Rsamtools", "extdata", "no_which_buffered_pileup.bam")

test_no_which_whole_file <- function() {
    bf <- BamFile(no_which_whole_file)
    xx <- pileup(bf)
    seqnames <- c(rep("chr1", 9), rep("chr2", 5))
    expected <- .tadf(seqnames=seqnames,
                      pos=c(1, 2, 3, 4, 4, 5, 5, 6, 7, seq(1,5)),
                      strand=rep("+",14),
                      nucleotide=c("A","A","A","A","C","A","C","C","C",
                        "A","A","A","C","C"),
                      count=c(1, 1, 2, rep(1,11)))
    ##print(xx);## str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}

test_no_which_buffered_yieldSize1_PREVENT_ZERO_ROW <- function() {
    bf <- open(BamFile(no_which_buf_pileup, yieldSize=1))
    ## return 1
    res <- pileup(bf)
    seqnames <- rep("chr1", 2L)
    expected <- .tadf(seqnames=seqnames,
                      pos=       c(  1,  2),
                      strand=    c("+","+"),
                      nucleotide=c("A","A"),
                      count=     c(  1,  1))
    ##print(res);##str(res)
    ##print(expected);str(expected)
    checkIdentical(expected, res)

    ## return 2
    res <- pileup(bf)
    seqnames <- rep("chr1", 4L)
    expected <- .tadf(seqnames=seqnames,
                      pos=       c(  3,  3,  4,  4),
                      strand=  rep("+",4L),
                      nucleotide=c("A","C","A","C"),
                      count=     c(  2,  1,  2,  1))
    ##print(res);##str(res)
    ## print(expected);##str(expected)
    checkIdentical(expected, res)

    ## return 3
    res <- pileup(bf)
    seqnames <- rep("chr1", 8L)
    expected <- .tadf(seqnames=seqnames,
                      pos=       c(  5,  5,  6,  6,  7,  7,  8,  9),
                      strand=  rep("+",8L),
                      nucleotide=c("A","C","A","C","A","C","A","A"),
                      count=     c(  3,  1,  2,  1,  2,  1,  1,  1))
    ##print(res);##str(res)
    ##print(expected);##str(expected)
    checkIdentical(expected, res)

    ## return 4 (detects EOI, flushes)
    res <- pileup(bf)
    seqnames <- rep("chr2", 5L)
    expected <- .tadf(seqnames=seqnames,
                      pos=         c(  5,  6,  7,  8,  9),
                      strand=    rep("+",5L),
                      nucleotide=rep("G",5L),
                      count=     rep(  1,5L))
    ##print(res);##str(res)
    ##print(expected);##str(expected)
    checkIdentical(expected, res)
}

## test_no_which_buffered_yieldSize1_NO_EAGER_PRINTING_DEMO <- function() {
##     bf <- open(BamFile(no_which_buf_pileup, yieldSize=1))
##     for(i in seq_len(5)) {
##         res <- pileup(bf)
##         message(paste0("printing return ", i))
##         print(res)
##         message()
##     }
##     res <- pileup(bf)
##     message(paste0("printing return 6"))
##     print(res);
## }
##test_no_which_buffered_yieldSize1_NO_EAGER_PRINTING_DEMO()

## Must be run with valgrind in order to verify
## test_no_which_buffered_finalizer <- function() {
##     bf <- open(BamFile(no_which_buf_pileup, yieldSize=1))
##     res <- pileup(bf)
##     bf <- open(BamFile(no_which_buf_pileup, yieldSize=1))
## }
##test_no_which_buffered_finalizer()
