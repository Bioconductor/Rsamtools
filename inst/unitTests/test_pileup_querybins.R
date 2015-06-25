suppressMessages({
    library(Rsamtools)
    library(RUnit)
})

## .run_all <- function() {
##     tests <- ls(envir=parent.frame())[grep('test_*', ls(envir=parent.frame()))]
##     cat(tests)
##     lapply(tests, function(x) { cat(paste0(x, "...\n")); do.call(x, list()) } )
## }

fl <- system.file(package="Rsamtools", "extdata", "querybins.bam")
bf <- BamFile(fl)

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
    c("A", "C", "G", "T", "N", "=", "-", "+")
}

.querybins.sam_seqlevels <- function() {
    fl <- system.file(package="Rsamtools", "extdata", "querybins.bam")
    bf <- BamFile(fl)
    seqlevels(bf)
}

## target as data.frame
.tadf <- function(pos=NULL, strand=NULL, nucleotide=NULL, left_bin=NULL,
                  query_bin=NULL,
                  count=NULL, which_label=NULL, seqnames=NULL,
                  use_ex1.bam_levels=FALSE) {
    if(is.null(pos) || is.null(count) || is.null(seqnames))
        stop("'pos', 'count', and 'seqnames' must not be 'NULL'")
    target <- data.frame(pos=as.integer(pos), stringsAsFactors=FALSE)
    seqnames_levels <- .querybins.sam_seqlevels()
    target <- cbind(seqnames=factor(seqnames, levels=seqnames_levels), target)

    if(!is.null(strand))
        target <- cbind(target, strand=factor(strand, levels=.s_levels()))
    if(!is.null(nucleotide))
        target <- cbind(target,
                        nucleotide=factor(nucleotide, levels=.n_levels()))
    if(!is.null(left_bin))
        target <- cbind(target, left_bin=left_bin)
    else if(!is.null(query_bin))
        target <- cbind(target, query_bin=query_bin)

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

##  1  2  3  4  5  6  7  8  9 10 11 12
##+ A  C  G  T  =  N
##-                   A  C  G  T  =  N
.sseq <- function() c("A", "C", "G", "T", "=", "N")
.s <- function()
    ScanBamParam(which=GRanges("sim", IRanges(1,15)))
.spos <- function() {
    ScanBamParam(which=GRanges("sim", IRanges(1,15)),
                 flag=scanBamFlag(isMinusStrand=FALSE))
}
.sneg <- function() {
    ScanBamParam(which=GRanges("sim", IRanges(1,15)),
                 flag=scanBamFlag(isMinusStrand=TRUE))
}

## distinguish_nucleotides=FALSE, ignore_query_Ns=FALSE
.pp <- function(qbins=NULL, lbins=NULL) {
    if(!is.null(qbins))
        PileupParam(distinguish_strands=FALSE, ignore_query_Ns=FALSE,
                    query_bins=qbins)
    else if(!is.null(lbins))
        PileupParam(distinguish_strands=FALSE, ignore_query_Ns=FALSE,
                    left_bins=lbins)
    else
        PileupParam(ignore_query_Ns=FALSE)
}

## PLUS STRAND TESTS

## query bin on plus strand, with finite (no Inf) positive bins
test_qbin_plus_pos_finite <- function() {
    sbp <- .spos()
    pp <- .pp(qbins=c(0, 3, 6))
    xx <- pileup(bf, scanBamParam=sbp, pileupParam=pp)
    seqnames <- space(bamWhich(sbp))
    expected <- .tadf(seqnames=seqnames,
                      pos=seq(1, 6),
                      nucleotide=.sseq(),
                      query_bin=c(rep("(0,3]", 3), rep("(3,6]", 3)),
                      count=rep(1, 6),
                      which_label=.mwls(sbp, length(seqnames)))
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_qbin_plus_pos_finite()

## query bin on plus strand, with finite (no Inf) positive bins, drop
## first and last cycle
test_qbin_plus_pos_finite_ends_dropped <- function() {
    sbp <- .spos()
    pp <- .pp(qbins=c(1, 3, 5))
    xx <- pileup(bf, scanBamParam=sbp, pileupParam=pp)
    seqnames <- space(bamWhich(sbp))
    expected <- .tadf(seqnames=seqnames,
                      pos=seq(2, 5),
                      nucleotide=c("C", "G", "T", "="),
                      query_bin=c("(1,3]", "(1,3]", "(3,5]", "(3,5]"),
                      count=rep(1, 4),
                      which_label=.mwls(sbp, length(seqnames)))
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_qbin_plus_pos_finite_ends_dropped()

## query bin on plus strand, with infinite positive bins
test_qbin_plus_pos_infinite <- function() {
    sbp <- .spos()
    pp <- .pp(qbins=c(0, 3, Inf))
    xx <- pileup(bf, scanBamParam=sbp, pileupParam=pp)
    seqnames <- space(bamWhich(sbp))
    expected <- .tadf(seqnames=seqnames,
                      pos=seq(1, 6),
                      nucleotide=.sseq(),
                      query_bin=c(rep("(0,3]", 3), rep("(3,Inf]", 3)),
                      count=rep(1,6),
                      which_label=.mwls(sbp, length(seqnames)))
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_qbin_plus_pos_infinite()

## query bin on plus strand, with finite (no -Inf) negative bins
test_qbin_plus_neg_finite <- function() {
    sbp <- .spos()
    pp <- .pp(qbins=c(-7, -4, -1))
    xx <- pileup(bf, scanBamParam=sbp, pileupParam=pp)
    seqnames <- space(bamWhich(sbp))
    expected <- .tadf(seqnames=seqnames,
                      pos=seq(1, 6),
                      nucleotide=.sseq(),
                      query_bin=c(rep("(-7,-4]", 3), rep("(-4,-1]", 3)),
                      count=rep(1,6),
                      which_label=.mwls(sbp, length(seqnames)))
    ## factor levels come back in different order, but values the same
    xx[["query_bin"]] <- as.character(xx[["query_bin"]])
    expected[["query_bin"]] <- as.character(expected[["query_bin"]])
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_qbin_plus_neg_finite()

## query bin on plus strand, with finite (no -Inf) negative bins, drop
## first and last cycle
test_qbin_plus_neg_finite_ends_dropped <- function() {
    sbp <- .spos()
    pp <- .pp(qbins=c(-6, -4, -2))
    xx <- pileup(bf, scanBamParam=sbp, pileupParam=pp)
    seqnames <- space(bamWhich(sbp))
    expected <- .tadf(seqnames=seqnames,
                      pos=seq(2, 5),
                      nucleotide=c("C", "G", "T", "="),
                      query_bin=c("(-6,-4]", "(-6,-4]", "(-4,-2]", "(-4,-2]"),
                      count=rep(1,4),
                      which_label=.mwls(sbp, length(seqnames)))
    ## factor levels come back in different order, but values the same
    xx[["query_bin"]] <- as.character(xx[["query_bin"]])
    expected[["query_bin"]] <- as.character(expected[["query_bin"]])
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_qbin_plus_neg_finite_ends_dropped()

## query bin on plus strand, with infinite negative bins
test_qbin_plus_neg_infinite <- function() {
    sbp <- .spos()
    pp <- .pp(qbins=c(-Inf, -4, -1))
    xx <- pileup(bf, scanBamParam=sbp, pileupParam=pp)
    seqnames <- space(bamWhich(sbp))
    expected <- .tadf(seqnames=seqnames,
                      pos=seq(1, 6),
                      nucleotide=.sseq(),
                      query_bin=c(rep("(-Inf,-4]", 3), rep("(-4,-1]", 3)),
                      count=rep(1,6),
                      which_label=.mwls(sbp, length(seqnames)))
    ## factor levels come back in different order, but values the same
    xx[["query_bin"]] <- as.character(xx[["query_bin"]])
    expected[["query_bin"]] <- as.character(expected[["query_bin"]])
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_qbin_plus_neg_infinite()

## MINUS STRAND TESTS

## query bin on minus strand, with finite (no Inf) positive bins
test_qbin_minus_pos_finite <- function() {
    sbp <- .sneg()
    pp <- .pp(qbins=c(0, 3, 6))
    xx <- pileup(bf, scanBamParam=sbp, pileupParam=pp)
    seqnames <- space(bamWhich(sbp))
    expected <- .tadf(seqnames=seqnames,
                      pos=seq(7, 12),
                      nucleotide=.sseq(),
                      query_bin=c(rep("(3,6]", 3), rep("(0,3]", 3)),
                      count=rep(1, 6),
                      which_label=.mwls(sbp, length(seqnames)))
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_qbin_minus_pos_finite()

## query bin on minus strand, with finite (no Inf) positive bins, drop
## first and last cycle
test_qbin_minus_pos_finite_ends_dropped <- function() {
    sbp <- .sneg()
    pp <- .pp(qbins=c(1, 3, 5))
    xx <- pileup(bf, scanBamParam=sbp, pileupParam=pp)
    seqnames <- space(bamWhich(sbp))
    expected <- .tadf(seqnames=seqnames,
                      pos=seq(8, 11),
                      nucleotide=c("C", "G", "T", "="),
                      query_bin=c("(3,5]", "(3,5]", "(1,3]", "(1,3]"),
                      count=rep(1, 4),
                      which_label=.mwls(sbp, length(seqnames)))
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_qbin_minus_pos_finite_ends_dropped()

## query bin on minus strand, with infinite positive bins
test_qbin_minus_pos_infinite <- function() {
    sbp <- .sneg()
    pp <- .pp(qbins=c(0, 3, Inf))
    xx <- pileup(bf, scanBamParam=sbp, pileupParam=pp)
    seqnames <- space(bamWhich(sbp))
    expected <- .tadf(seqnames=seqnames,
                      pos=seq(7, 12),
                      nucleotide=.sseq(),
                      query_bin=c(rep("(3,Inf]", 3), rep("(0,3]", 3)),
                      count=rep(1,6),
                      which_label=.mwls(sbp, length(seqnames)))
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_qbin_minus_pos_infinite()

## query bin on minus strand, with finite (no -Inf) negative bins
test_qbin_minus_neg_finite <- function() {
    sbp <- .sneg()
    pp <- .pp(qbins=c(-7, -4, -1))
    xx <- pileup(bf, scanBamParam=sbp, pileupParam=pp)
    seqnames <- space(bamWhich(sbp))
    expected <- .tadf(seqnames=seqnames,
                      pos=seq(7, 12),
                      nucleotide=.sseq(),
                      query_bin=c(rep("(-4,-1]", 3), rep("(-7,-4]", 3)),
                      count=rep(1,6),
                      which_label=.mwls(sbp, length(seqnames)))
    ## factor levels come back in different order, but values the same
    xx[["query_bin"]] <- as.character(xx[["query_bin"]])
    expected[["query_bin"]] <- as.character(expected[["query_bin"]])
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_qbin_minus_neg_finite()

## query bin on minus strand, with finite (no -Inf) negative bins, drop
## first and last cycle
test_qbin_minus_neg_finite_ends_dropped <- function() {
    sbp <- .sneg()
    pp <- .pp(qbins=c(-6, -4, -2))
    xx <- pileup(bf, scanBamParam=sbp, pileupParam=pp)
    seqnames <- space(bamWhich(sbp))
    expected <- .tadf(seqnames=seqnames,
                      pos=seq(8, 11),
                      nucleotide=c("C", "G", "T", "="),
                      query_bin=c("(-4,-2]", "(-4,-2]", "(-6,-4]", "(-6,-4]"),
                      count=rep(1,4),
                      which_label=.mwls(sbp, length(seqnames)))
    ## factor levels come back in different order, but values the same
    xx[["query_bin"]] <- as.character(xx[["query_bin"]])
    expected[["query_bin"]] <- as.character(expected[["query_bin"]])
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_qbin_minus_neg_finite_ends_dropped()

## query bin on minus strand, with infinite negative bins
test_qbin_minus_neg_infinite <- function() {
    sbp <- .sneg()
    pp <- .pp(qbins=c(-Inf, -4, -1))
    xx <- pileup(bf, scanBamParam=sbp, pileupParam=pp)
    seqnames <- space(bamWhich(sbp))
    expected <- .tadf(seqnames=seqnames,
                      pos=seq(7, 12),
                      nucleotide=.sseq(),
                      query_bin=c(rep("(-4,-1]", 3), rep("(-Inf,-4]", 3)),
                      count=rep(1,6),
                      which_label=.mwls(sbp, length(seqnames)))
    ## factor levels come back in different order, but values the same
    xx[["query_bin"]] <- as.character(xx[["query_bin"]])
    expected[["query_bin"]] <- as.character(expected[["query_bin"]])
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_qbin_minus_neg_infinite()

## COMBINED STRANDS

## query bin on both strands, pos bins, with 3 bins
test_qbin_both_pos_3bin <- function() {
    sbp <- .s()
    pp <- .pp(qbins=c(0, 2, 4, 6))
    xx <- pileup(bf, scanBamParam=sbp, pileupParam=pp)
    seqnames <- space(bamWhich(sbp))
    plusbins <- c("(0,2]","(0,2]","(2,4]","(2,4]","(4,6]","(4,6]")
    expected <- .tadf(seqnames=seqnames,
                      pos=seq(1, 12),
                      nucleotide=c(.sseq(), .sseq()),
                      query_bin=c(plusbins, rev(plusbins)),
                      count=rep(1,12),
                      which_label=.mwls(sbp, length(seqnames)))
    ## factor levels come back in different order, but values the same
    xx[["query_bin"]] <- as.character(xx[["query_bin"]])
    expected[["query_bin"]] <- as.character(expected[["query_bin"]])
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_qbin_both_pos_3bin()

## query bin on both strands, neg bins, with 3 bins
test_qbin_both_neg_3bin <- function() {
    sbp <- .s()
    pp <- .pp(qbins=c(-7, -5, -3, -1))
    xx <- pileup(bf, scanBamParam=sbp, pileupParam=pp)
    seqnames <- space(bamWhich(sbp))
    plusbins <- c("(-7,-5]","(-7,-5]","(-5,-3]","(-5,-3]","(-3,-1]","(-3,-1]")
    expected <- .tadf(seqnames=seqnames,
                      pos=seq(1, 12),
                      nucleotide=c(.sseq(), .sseq()),
                      query_bin=c(plusbins, rev(plusbins)),
                      count=rep(1,12),
                      which_label=.mwls(sbp, length(seqnames)))
    ## factor levels come back in different order, but values the same
    xx[["query_bin"]] <- as.character(xx[["query_bin"]])
    expected[["query_bin"]] <- as.character(expected[["query_bin"]])
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_qbin_both_neg_3bin()

