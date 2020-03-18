suppressMessages({
    library(Rsamtools)
    library(RUnit)
})

## .run_all <- function() {
##     tests <- ls(envir=parent.frame())[grep('test_*', ls(envir=parent.frame()))]
##     cat(tests)
##     lapply(tests, function(x) { cat(paste0(x, "...\n")); do.call(x, list()) } )
## }

fl <- system.file(package="Rsamtools", "extdata", "revbins.bam")
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

.revbins.sam_seqlevels <- function() {
    fl <- system.file(package="Rsamtools", "extdata", "revbins.bam")
    bf <- BamFile(fl)
    seqlevels(bf)
}

## target as data.frame
.tadf <- function(pos=NULL, strand=NULL, nucleotide=NULL, left_bin=NULL,
                  count=NULL, which_label=NULL, seqnames=NULL,
                  use_ex1.bam_levels=FALSE) {
    if(is.null(pos) || is.null(count) || is.null(seqnames))
        stop("'pos', 'count', and 'seqnames' must not be 'NULL'")
    target <- data.frame(pos=as.integer(pos))
    seqnames_levels <- .revbins.sam_seqlevels()
    target <- cbind(seqnames=factor(seqnames, levels=seqnames_levels), target)

    if(!is.null(strand))
        target <- cbind(target,
                        strand=factor(strand, levels=.s_levels()))
    if(!is.null(nucleotide))
        target <- cbind(target,
                        nucleotide=factor(nucleotide, levels=.n_levels()))
    if(!is.null(left_bin))
        target <- cbind(target,
                        left_bin=left_bin,
                        stringsAsFactors=TRUE)
    target <- cbind(target, count=as.integer(count))
    if(!is.null(which_label))
        target <- cbind(target, which_label=which_label, stringsAsFactors=TRUE)
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

.revbinsam_empty_df <- function(param) {
    wl <- .mwls(param, 0)
    .tadf(seqnames=character(),
          pos=integer(),
          strand=character(),
          nucleotide=character(),
          count=integer(),
          which_label=wl
          )
}

## 1-width
.simple <- function()
    ScanBamParam(which=GRanges("simple", IRanges(1,10)))
## 2-width staggered
.bistag <- function()
    ScanBamParam(which=GRanges("bistag", IRanges(1,10)))
## 3-width staggered
.tristag <- function()
    ScanBamParam(which=GRanges("tristag", IRanges(1,10)))
## for testing -Inf
.inf <- function()
    ScanBamParam(which=GRanges("inf", IRanges(1,10)))
## different strands
.diffstr <- function()
    ScanBamParam(which=GRanges("diffstr", IRanges(1,10)))

## single bin, width 1
test_revbin_last <- function() {
    sb_param <- .simple()
    pp <- PileupParam(
        distinguish_strands=FALSE,
        left_bins=c(-2, -1))
    xx <- pileup(bf, scanBamParam=sb_param, pileupParam=pp)
    seqnames <- space(bamWhich(sb_param))
    expected <- .tadf(seqnames=seqnames,
                      pos=seq(1L, 5L),
                      nucleotide=rep("=", 5),
                      left_bin=rep("(-2,-1]", 5),
                      count=rep(1L,5L),
                      which_label=.mwls(sb_param, length(seqnames)))
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_revbin_last()

## single bin, width 1
test_revbin_1_from_last <- function() {
    sb_param <- .simple()
    pp <- PileupParam(
        distinguish_strands=FALSE,
        left_bins=c(-3, -2))
    xx <- pileup(bf, scanBamParam=sb_param, pileupParam=pp)
    seqnames <- space(bamWhich(sb_param))
    expected <- .tadf(seqnames=seqnames,
                      pos=seq(1L, 4L),
                      nucleotide=rep("A", 4),
                      left_bin=rep("(-3,-2]", 4),
                      count=rep(1L,4L),
                      which_label=.mwls(sb_param, length(seqnames)))
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_revbin_1_from_last()

## single bin, width 1
test_revbin_2_from_last <- function() {
    sb_param <- .simple()
    pp <- PileupParam(
        distinguish_strands=FALSE,
        left_bins=c(-4, -3))
    xx <- pileup(bf, scanBamParam=sb_param, pileupParam=pp)
    seqnames <- space(bamWhich(sb_param))
    expected <- .tadf(seqnames=seqnames,
                      pos=seq(1L, 3L),
                      nucleotide=rep("C", 3),
                      left_bin=rep("(-4,-3]", 3),
                      count=rep(1L,3L),
                      which_label=.mwls(sb_param, length(seqnames)))
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_revbin_2_from_last()

## single bin, width 2
test_revbin_width_2 <- function() {
    sb_param <- .simple()
    pp <- PileupParam(
        distinguish_strands=FALSE,
        left_bins=c(-3, -1))
    xx <- pileup(bf, scanBamParam=sb_param, pileupParam=pp)
    seqnames <- space(bamWhich(sb_param))
    expected <- .tadf(seqnames=seqnames,
                      pos=c(1, 1, 2, 2, 3, 3, 4, 4, 5),
                      nucleotide=c("=", "A", "=", "A", "=", "A", "=", "A", "="),
                      left_bin=rep("(-3,-1]", 9),
                      count=rep(1L,9L),
                      which_label=.mwls(sb_param, length(seqnames)))
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_revbin_width_2()

## 2 bin, width 1, all satisfy length
test_revbin_2_bin_width_1_321 <- function() {
    sb_param <- .simple()
    pp <- PileupParam(
        distinguish_strands=FALSE,
        left_bins=c(-3, -2, -1))
    xx <- pileup(bf, scanBamParam=sb_param, pileupParam=pp)
    seqnames <- space(bamWhich(sb_param))
    expected <- .tadf(seqnames=seqnames,
                      pos=c(1, 1, 2, 2, 3, 3, 4, 4, 5),
                      nucleotide=c("=", "A", "=", "A", "=", "A", "=", "A", "="),
                      left_bin=c("(-2,-1]","(-3,-2]","(-2,-1]","(-3,-2]",
                                  "(-2,-1]","(-3,-2]","(-2,-1]","(-3,-2]",
                                  "(-2,-1]"),
                      count=rep(1L,9L),
                      which_label=.mwls(sb_param, length(seqnames)))
    ## factor levels come back in different order, but values the same
    xx[["left_bin"]] <- as.character(xx[["left_bin"]])
    expected[["left_bin"]] <- as.character(expected[["left_bin"]])
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_revbin_2_bin_width_1_321()

## 2 bin, width 1, all but one satisfy length
test_revbin_2_bin_width_1_432 <- function() {
    sb_param <- .simple()
    pp <- PileupParam(
        distinguish_strands=FALSE,
        left_bins=c(-4, -3, -2))
    xx <- pileup(bf, scanBamParam=sb_param, pileupParam=pp)
    seqnames <- space(bamWhich(sb_param))
    expected <- .tadf(seqnames=seqnames,
                      pos=c(1, 1, 2, 2, 3, 3, 4),
                      nucleotide=c("A", "C", "A", "C", "A", "C", "A"),
                      left_bin=c("(-3,-2]","(-4,-3]","(-3,-2]", "(-4,-3]",
                                  "(-3,-2]","(-4,-3]","(-3,-2]"),
                      count=rep(1L,7L),
                      which_label=.mwls(sb_param, length(seqnames)))
    ## factor levels come back in different order, but values the same
    xx[["left_bin"]] <- as.character(xx[["left_bin"]])
    expected[["left_bin"]] <- as.character(expected[["left_bin"]])
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_revbin_2_bin_width_1_432()

## arbitrary bins
test_revbin_arb_bins <- function() {
    sb_param <- .simple()
    pp <- PileupParam(
        distinguish_strands=FALSE,
        left_bins=c(-5, -4, -3, -1))
    xx <- pileup(bf, scanBamParam=sb_param, pileupParam=pp)
    seqnames <- space(bamWhich(sb_param))
    expected <-
        .tadf(seqnames=seqnames,
            pos=c(rep(1, 4), rep(2, 4), rep(3, 3), 4, 4, 5),
            nucleotide=c("=","A","C","G","=","A","C","G","=","A","C","=","A","="),
            left_bin=c("(-3,-1]","(-3,-1]","(-4,-3]","(-5,-4]",
                        "(-3,-1]","(-3,-1]","(-4,-3]","(-5,-4]",
                        "(-3,-1]","(-3,-1]","(-4,-3]",
                        "(-3,-1]","(-3,-1]","(-3,-1]"),
            count=rep(1L,14L),
            which_label=.mwls(sb_param, length(seqnames)))
    ## factor levels come back in different order, but values the same
    xx[["left_bin"]] <- as.character(xx[["left_bin"]])
    expected[["left_bin"]] <- as.character(expected[["left_bin"]])
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_revbin_arb_bins()

## arbitrary bins, collapse nucleotides
test_revbin_arb_bins_collapse_nucs <- function() {
    sb_param <- .simple()
    pp <- PileupParam(
        distinguish_strands=FALSE,
        distinguish_nucleotides=FALSE,
        left_bins=c(-5, -4, -3, -1))
    xx <- pileup(bf, scanBamParam=sb_param, pileupParam=pp)
    seqnames <- space(bamWhich(sb_param))
    expected <-
        .tadf(seqnames=seqnames,
            pos=c(rep(1, 3), rep(2, 3), rep(3, 2), 4, 5),
            left_bin=c("(-5,-4]","(-4,-3]","(-3,-1]",
                        "(-5,-4]","(-4,-3]","(-3,-1]",
                        "(-4,-3]","(-3,-1]",
                        "(-3,-1]",
                        "(-3,-1]"),
            count=c(1,1,2,1,1,2,1,2,2,1),
            which_label=.mwls(sb_param, length(seqnames)))
    ## factor levels come back in different order, but values the same
    xx[["left_bin"]] <- as.character(xx[["left_bin"]])
    expected[["left_bin"]] <- as.character(expected[["left_bin"]])
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_revbin_arb_bins_collapse_nucs()

## bin range to beginning, collapse nucleotides
test_revbin_collapse_to_beginning <- function() {
    sb_param <- .simple()
    pp <- PileupParam(
        distinguish_strands=FALSE,
        distinguish_nucleotides=FALSE,
        left_bins=c(-6, -2))
    xx <- pileup(bf, scanBamParam=sb_param, pileupParam=pp)
    seqnames <- space(bamWhich(sb_param))
    expected <-
        .tadf(seqnames=seqnames,
            pos=c(1, 2, 3, 4),
            left_bin=c(rep("(-6,-2]", 4)),
            count=c(4, 3, 2, 1),
            which_label=.mwls(sb_param, length(seqnames)))
    ## factor levels come back in different order, but values the same
    xx[["left_bin"]] <- as.character(xx[["left_bin"]])
    expected[["left_bin"]] <- as.character(expected[["left_bin"]])
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_revbin_collapse_to_beginning()

## staggered, reads length 2, separate bins end to beginning
test_revbin_staggered_sep_bins <- function() {
    sb_param <- .bistag()
    pp <- PileupParam(
        distinguish_strands=FALSE,
        ignore_query_Ns=FALSE,
        left_bins=c(-3, -2, -1))
    xx <- pileup(bf, scanBamParam=sb_param, pileupParam=pp)
    seqnames <- space(bamWhich(sb_param))
    expected <-
        .tadf(seqnames=seqnames,
            pos=c(1, 2, 2, 3, 3, 4, 4, 5, 5, 6),
            nucleotide=c("A","=","C","=","G","=","T","=","N","="),
            left_bin=rep(c("(-3,-2]", "(-2,-1]"), 5),
            count=rep(1, 10),
            which_label=.mwls(sb_param, length(seqnames)))
    ## factor levels come back in different order, but values the same
    xx[["left_bin"]] <- as.character(xx[["left_bin"]])
    expected[["left_bin"]] <- as.character(expected[["left_bin"]])
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_revbin_staggered_sep_bins()

## staggered, reads length 2, single bin end to beginning
test_revbin_staggered_single_bin <- function() {
    sb_param <- .bistag()
    pp <- PileupParam(
        distinguish_strands=FALSE,
        ignore_query_Ns=FALSE,
        left_bins=c(-3, -1))
    xx <- pileup(bf, scanBamParam=sb_param, pileupParam=pp)
    seqnames <- space(bamWhich(sb_param))
    expected <-
        .tadf(seqnames=seqnames,
            pos=c(1, 2, 2, 3, 3, 4, 4, 5, 5, 6),
            nucleotide=c("A","=","C","=","G","=","T","=","N","="),
            left_bin=rep("(-3,-1]", 10),
            count=rep(1, 10),
            which_label=.mwls(sb_param, length(seqnames)))
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_revbin_staggered_single_bin()

## staggered, reads length 2, single bin end to beginning, collapse nucs
test_revbin_staggered_single_bin_collapse_nucs <- function() {
    sb_param <- .bistag()
    pp <- PileupParam(
        distinguish_strands=FALSE,
        distinguish_nucleotides=FALSE,
        ignore_query_Ns=FALSE,
        left_bins=c(-3, -1))
    xx <- pileup(bf, scanBamParam=sb_param, pileupParam=pp)
    seqnames <- space(bamWhich(sb_param))
    expected <-
        .tadf(seqnames=seqnames,
            pos=c(1, 2, 3, 4, 5, 6),
            left_bin=rep("(-3,-1]", 6),
            count=c(1, 2, 2, 2, 2, 1),
            which_label=.mwls(sb_param, length(seqnames)))
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_revbin_staggered_single_bin_collapse_nucs()

## staggered, reads length 2-3, c(-4, -2)
test_revbin_staggered_omit_end <- function() {
    sb_param <- .tristag()
    pp <- PileupParam(
        distinguish_strands=FALSE,
        ignore_query_Ns=FALSE,
        left_bins=c(-4, -2))
    xx <- pileup(bf, scanBamParam=sb_param, pileupParam=pp)
    seqnames <- space(bamWhich(sb_param))
    expected <-
        .tadf(seqnames=seqnames,
            pos=c(1, 2, 3, 3, 4, 4, 5, 5, 6),
            nucleotide=c("A","C","=","G","=","T","=","N","="),
            left_bin=rep("(-4,-2]", 9),
            count=rep(1, 9),
            which_label=.mwls(sb_param, length(seqnames)))
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_revbin_staggered_omit_end()

## staggered, reads length 2-3, c(-4, -2)
test_revbin_staggered_omit_end_collapse_nucs <- function() {
    sb_param <- .tristag()
    pp <- PileupParam(
        distinguish_strands=FALSE,
        distinguish_nucleotides=FALSE,
        ignore_query_Ns=FALSE,
        left_bins=c(-4, -2))
    xx <- pileup(bf, scanBamParam=sb_param, pileupParam=pp)
    seqnames <- space(bamWhich(sb_param))
    expected <-
        .tadf(seqnames=seqnames,
            pos=c(1, 2, 3, 4, 5, 6),
            left_bin=rep("(-4,-2]", 6),
            count=c(1, 1, 2, 2, 2, 1),
            which_label=.mwls(sb_param, length(seqnames)))
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_revbin_staggered_omit_end_collapse_nucs()

## -INF

## c(-Inf, -3, -2)
test_revbin_inf <- function() {
    sb_param <- .inf()
    pp <- PileupParam(
        distinguish_strands=FALSE,
        left_bins=c(-Inf, -3, -2))
    xx <- pileup(bf, scanBamParam=sb_param, pileupParam=pp)
    seqnames <- space(bamWhich(sb_param))
    expected <-
        .tadf(seqnames=seqnames,
            pos=c(1, 2, 3, 3, 4, 4, 5, 5, 6),
            nucleotide=c("A", "A", "A", "C", "A", "C", "A", "C", "C"),
            left_bin=c(rep("(-Inf,-3]", 6), "(-3,-2]", "(-Inf,-3]", "(-3,-2]"),
            count=rep(1, 9),
            which_label=.mwls(sb_param, length(seqnames)))
    ## factor levels come back in different order, but values the same
    xx[["left_bin"]] <- as.character(xx[["left_bin"]])
    expected[["left_bin"]] <- as.character(expected[["left_bin"]])
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_revbin_inf()

## c(-Inf, -3, -2), collapse nucs
test_revbin_inf_collapse_nucs <- function() {
    sb_param <- .inf()
    pp <- PileupParam(
        distinguish_strands=FALSE,
        distinguish_nucleotides=FALSE,
        left_bins=c(-Inf, -3, -2))
    xx <- pileup(bf, scanBamParam=sb_param, pileupParam=pp)
    seqnames <- space(bamWhich(sb_param))
    expected <-
        .tadf(seqnames=seqnames,
            pos=c(1, 2, 3, 4, 5, 5, 6),
            left_bin=c(rep("(-Inf,-3]", 5), "(-3,-2]", "(-3,-2]"),
            count=c(1, 1, 2, 2, 1, 1, 1),
            which_label=.mwls(sb_param, length(seqnames)))
    ## factor levels come back in different order, but values the same
    xx[["left_bin"]] <- as.character(xx[["left_bin"]])
    expected[["left_bin"]] <- as.character(expected[["left_bin"]])
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_revbin_inf_collapse_nucs()

## end -Inf

## DIFFERENT STRANDS

## c(-4, -1), dist strands TRUE, dist nucs TRUE
test_revbin_diff_strands_distinguish_strands <- function() {
    sb_param <- .diffstr()
    pp <- PileupParam(
        ignore_query_Ns=FALSE,
        left_bins=c(-4, -1))
    xx <- pileup(bf, scanBamParam=sb_param, pileupParam=pp)
    seqnames <- space(bamWhich(sb_param))
    expected <-
        .tadf(seqnames=seqnames,
            pos=c(1, 2, 2, 3, 3, 4),
            strand=c("+", "+", "-", "-", "+", "-"),
            nucleotide=c("A", "C", "T", "=", "G", "N"),
            left_bin=rep("(-4,-1]", 6),
            count=rep(1, 6),
            which_label=.mwls(sb_param, length(seqnames)))
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_revbin_diff_strands_distinguish_strands()

## c(-4, -1), dist strands FALSE, dist nucs FALSE
test_revbin_diff_strands_collapse_all <- function() {
    sb_param <- .diffstr()
    pp <- PileupParam(
        distinguish_strands=FALSE,
        distinguish_nucleotides=FALSE,
        ignore_query_Ns=FALSE,
        left_bins=c(-4, -1))
    xx <- pileup(bf, scanBamParam=sb_param, pileupParam=pp)
    seqnames <- space(bamWhich(sb_param))
    expected <-
        .tadf(seqnames=seqnames,
            pos=c(1, 2, 3, 4),
            left_bin=rep("(-4,-1]", 4),
            count=c(1, 2, 2, 1),
            which_label=.mwls(sb_param, length(seqnames)))
    ##print(xx); ##str(xx);
    ##print(expected); ##str(expected)
    checkIdentical(expected, xx)
}
##test_revbin_diff_strands_collapse_all()
