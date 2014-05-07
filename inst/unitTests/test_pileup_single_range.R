suppressMessages({
    library(Rsamtools)
    library(RUnit)
})

## .run_all <- function() {
##     tests <- ls(envir=parent.frame())[grep('test_*', ls(envir=parent.frame()))]
##     cat(tests)
##     lapply(tests, function(x) { cat(paste0(x, "...\n")); do.call(x, list()) } )
## }

fl <- system.file(package="Rsamtools", "extdata", "tiny.bam")
bf <- BamFile(fl)

.unordered_check <- function(expected, xx) {
    ##xx <- as.data.frame(xx) ## no longer need coercion
    xx <- xx[do.call(order,xx),]
    rownames(xx) <- NULL
    ## expected <- as.data.frame(expected) ## no longer need coercion
    expected <- expected[do.call(order,expected),]
    rownames(expected) <- NULL
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

## target as data.frame
.tadf <- function(pos=NULL, strand=NULL, nucleotide=NULL, count=NULL,
                  which_label=NULL, seqnames=NULL, use_ex1.bam_levels=FALSE) {
    if(is.null(pos) || is.null(count) || is.null(seqnames))
        stop("'pos', 'count', and 'seqnames' must not be 'NULL'")
    target <- data.frame(pos=as.integer(pos), stringsAsFactors=FALSE)
    seqnames_levels <- .tiny.sam_seqlevels()
    if(use_ex1.bam_levels)
        seqnames_levels <- .ex1.sam_seqlevels()        
    target <- cbind(seqnames=factor(seqnames, levels=seqnames_levels), target)
    
    if(!is.null(strand))
        target <- cbind(target, strand=factor(strand, levels=.s_levels()))
    if(!is.null(nucleotide))
        target <- cbind(target,
                        nucleotide=factor(nucleotide, levels=.n_levels()))
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

.tinysam_empty_df <- function(param) {
    wl <- .mwls(param, 0)
    .tadf(seqnames=character(),
          pos=integer(),
          strand=character(),
          nucleotide=character(),
          count=integer(),
          which_label=wl
          )
}

.seq1 <- function()
    ScanBamParam(which=GRanges("seq1", IRanges(1,10)))
.seq2 <- function()
    ScanBamParam(which=GRanges("seq2", IRanges(1,10)))
.diff_strands <- function()
    ScanBamParam(which=GRanges("diff_strands", IRanges(1,10)))
.coll_nucs_basic <- function()
    ScanBamParam(which=GRanges("coll_nucs_basic", IRanges(1,10)))
.coll_nucs_adv <- function()
    ScanBamParam(which=GRanges("coll_nucs_adv", IRanges(1,10)))
.dist_all <- function()
    ScanBamParam(which=GRanges("dist_all", IRanges(1,10)))
.min_mapq <- function()
    ScanBamParam(which=GRanges("min_mapq", IRanges(1,10)))
.min_base_qual <- function()
    ScanBamParam(which=GRanges("min_base_qual", IRanges(1,10)))
.max_depth <- function()
    ScanBamParam(which=GRanges("max_depth", IRanges(1,10)))
.full_seq_range <- function()
    ScanBamParam(which=GRanges("full_seq_range", IRanges(1,28)))
.Biostrings_DNA_ALPHABET <- function()
    ScanBamParam(which=GRanges("Biostrings_DNA_ALPHABET", IRanges(1,18)))
.ext_cigar_hard <- function()
    ScanBamParam(which=GRanges("ext_cigar_hard", IRanges(1,10)))
.ext_cigar_soft <- function()
    ScanBamParam(which=GRanges("ext_cigar_soft", IRanges(1,10)))
.ref_skip <- function()
    ScanBamParam(which=GRanges("ref_skip", IRanges(1,10)))
.include_deletions <- function()
    ScanBamParam(which=GRanges("include_deletions", IRanges(1,10)))
.min_nuc_depth_mono <- function()
    ScanBamParam(which=GRanges("min_nuc_depth_mono", IRanges(1,10)))
.min_nuc_depth_poly <- function()
    ScanBamParam(which=GRanges("min_nuc_depth_poly", IRanges(1,10)))
.min_minor_allele_depth <- function()
    ScanBamParam(which=GRanges("min_minor_allele_depth", IRanges(1,10)))
.multi_range_single_rname <- function()
    ScanBamParam(which=GRanges("multi_range_single_rname", IRanges(c(1,3),width=2)))
.filter_strands <- function()
    ScanBamParam(which=GRanges("filter_strands", IRanges(1, 10)))
.seqnames_from_tinysam_c <- function()
    ScanBamParam(which=GRanges("seqnames_from_tinysam_c", IRanges(1,2)))
.all_nuc_levels <- function()
    ScanBamParam(which=GRanges("all_nuc_levels", IRanges(1, 7)))

## test all nuc levels represented
test_all_nuc_levels <- function() {
    scanBamParam <- .all_nuc_levels()
    pileupParam <- PileupParam(
        distinguish_strands=FALSE,
        ignore_query_Ns=FALSE,
        include_deletions=TRUE)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam))
    expected <- .tadf(seqnames=seqnames,
                      pos=seq(1, 7),
                      count=rep(1, 7),
                      nucleotide=c("A", "C", "-", "G", "T", "N", "="),
                      which_label=.mwls(scanBamParam, length(seqnames)))
    ##print(xx);str(xx);print(expected);str(expected)
    checkIdentical(expected, xx)
}

## test that seqnames are set from C correctly
test_seqnames_from_tinysam_in_c <- function() {
    scanBamParam <- .seqnames_from_tinysam_c()
    pileupParam <- PileupParam(
                               distinguish_strands=FALSE,
                               distinguish_nucleotides=FALSE)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam))
    expected <- .tadf(seqnames=seqnames,
                      pos = 1,
                      count = 1,
                      which_label=.mwls(scanBamParam, length(seqnames)))
    ##print(xx);str(xx);print(expected);str(expected)
    checkIdentical(expected, xx)
}

test_scanBam_filter_criteria_dropsome <- function() {
    ## assume that if one filtering criterion works, all work
    scanBamParam <- .filter_strands()
    flag <- scanBamFlag(isMinusStrand=TRUE)
    bamFlag(scanBamParam) <- flag
    pileupParam <- PileupParam(distinguish_nucleotides=FALSE)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam))
    expected <- .tadf(seqnames=seqnames,
                      pos   = 1,
                      strand= "-",
                      count = 1,
                      which_label=.mwls(scanBamParam, length(seqnames)))
    ## print(xx);str(xx)
    ## print(expected);str(expected)
    checkIdentical(expected, xx)
}

test_scanBam_filter_criteria_dropnone <- function() {
    ## assume that if one filtering criterion works, all work
    scanBamParam <- .filter_strands()
    pileupParam <- PileupParam(distinguish_nucleotides=FALSE)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam))
    expected <- .tadf(seqnames=seqnames,
                      pos = c(1, 1),
                      strand = c("+", "-"),
                      count = c(1, 1),
                      which_label=.mwls(scanBamParam, length(seqnames)))
    ## print(xx);str(xx)
    ## print(expected);str(expected)
    checkIdentical(expected, xx)
}

test_data.frame_format <- function() {
    fl <- system.file(package="Rsamtools", "extdata", "ex1.bam")
    bf <- BamFile(fl)
    scanBamParam <- ScanBamParam(which=GRanges(c("seq1", "seq2"), IRanges(1, 3)))
    pileupParam <- PileupParam(distinguish_strands=TRUE,
                               distinguish_nucleotides=TRUE)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    expected <- .tadf(
        seqnames=c(rep("seq1", 3), rep("seq2", 6)),
        pos=c(1, 2, 3, 1, 1, 2, 2, 3, 3),
        strand=c("+", "+", "+", "+", "-", "+", "-", "+", "-"),
        nucleotide=c("C", "A", "C", "T", "T", "T", "T", "C", "C"),
        count=c(1, 1, 2, 2, 1, 3, 1, 3, 1),
        which_label=c(rep("seq1:1-3", 3), rep("seq2:1-3", 6)),
        use_ex1.bam_levels=TRUE)
    ##print(xx);str(xx);print(expected);str(expected)
    checkIdentical(expected, xx)
}

## FIX ME: better name
test_multirange_yield_clear <- function() {
    scanBamParam <- .multi_range_single_rname()
    pileupParam <- PileupParam(distinguish_strands=FALSE,
                               distinguish_nucleotides=TRUE)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam))
    expected <- .tadf(seqnames=seqnames,
                      pos = seq(1L, 4L),
                      nucleotide=c("A", "C", "G", "T"),
                      count=rep(1L, 4L),
                      which_label=.mwls(scanBamParam, c(2, 2)))
    ##print(xx);str(xx);print(expected);str(expected)
    checkIdentical(expected, xx)
}

## MULTI-READ DISTINGUISH*

test_distinguish_none <- function() {
    scanBamParam <- .dist_all()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=FALSE)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos=     1L,
                     count=   20L,
                     which_label=.mwls(scanBamParam, length(seqnames)))
    ##print(xx);print(expected)
    checkIdentical(expected, xx)
}

test_distinguish_nuc_only <- function() {
    scanBamParam <- .dist_all()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=TRUE)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos=     rep(1L, 4L),
                     nucleotide=c("A", "C", "G", "T"),
                     count=   rep(5L, 4L),
                     which_label=.mwls(scanBamParam, length(seqnames)))
    ##print(xx);print(expected)
    .unordered_check(expected, xx)
}
test_distinguish_strand_only <- function() {
    scanBamParam <- .dist_all()
    pileupParam <- PileupParam(
                              distinguish_strands=TRUE,
                              distinguish_nucleotides=FALSE)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos=   c(1L, 1L),
                     strand=c("+", "-"),
                     count= c(10L, 10L),
                      which_label=.mwls(scanBamParam, length(seqnames)))
    ##print(xx);print(expected)
    checkIdentical(expected, xx)
}
test_distinguish_all <- function() {
    ## For all the marbles
    scanBamParam <- .dist_all()
    pileupParam <- PileupParam(
                              distinguish_strands=TRUE,
                              distinguish_nucleotides=TRUE,
                              ignore_query_Ns=FALSE)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos=      rep( 1L,  8L),
                     strand=     c("+", "-", "+", "-", "+", "-", "+", "-"),
                     nucleotide= c("A", "C", "G", "T", "C", "G", "T", "A"),
                     count=      c( 3L,  3L,  3L,  3L,  2L,  2L,  2L,  2L),
                      which_label=.mwls(scanBamParam, length(seqnames)))
    .unordered_check(expected, xx)
}

## SINGLE-READ DISTINGUISH*

test_distinguish_nucleotides_adv <- function() {
    ## Test that a position with more than two alignments with diff nucleotides
    ## gets counted correctly
    scanBamParam <- .coll_nucs_adv()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=TRUE,
                              ignore_query_Ns=FALSE)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos=        c( 1L,  1L, 1L),
                     nucleotide= c("G", "A", "T"),
                     count=      c( 2L,  1L, 1L),
                      which_label=.mwls(scanBamParam, length(seqnames)))
    ##print(xx);print(expected)
    .unordered_check(expected, xx)
}

test_distinguish_nucleotides_basic <- function() {
    ## Test that a position with more than two alignments with *same* nucleotide
    ## gets counted correctly
    scanBamParam <- .coll_nucs_basic()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=TRUE,
                              ignore_query_Ns=TRUE)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos=        seq(1L, 5L, 1L),
                     nucleotide= rep("C", 5L),
                     count=      rep(2L, 5L),
                      which_label=.mwls(scanBamParam, length(seqnames)))
    checkIdentical(expected, xx)
}

test_no_distinguish_strands_mixed <- function() {
    ## NOTE: even though this looks like "distinguish nothing", it's not:
    ## test seqs in tiny.sam are identical POS and SEQ, but different strands
    scanBamParam <- .diff_strands()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=FALSE,
                              ignore_query_Ns=TRUE)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos=   seq(1L, 5L, 1L),
                     count= rep(2L, 5L),
                      which_label=.mwls(scanBamParam, length(seqnames)))
    checkIdentical(expected, xx)
}

test_distinguish_strands_mixed <- function() {
    scanBamParam <- .diff_strands()
    pileupParam <- PileupParam(
                              distinguish_strands=TRUE,
                              distinguish_nucleotides=FALSE,
                              ignore_query_Ns=FALSE)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos=   c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 5L),
                     strand=rep(c("+", "-"), 5),
                     count= rep(1L, 10),
                      which_label=.mwls(scanBamParam, length(seqnames)))
    checkIdentical(expected, xx)
}

## MIN MINOR ALLELE DEPTH

test_minor_allele_depth_dropnone <- function() {
    scanBamParam <- .min_minor_allele_depth()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=FALSE,
                              min_nucleotide_depth=0L,
                              min_minor_allele_depth=1L)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos=  1L,
                     count=4L,
                      which_label=.mwls(scanBamParam, length(seqnames)))
    ##print(xx);print(expected)
    checkIdentical(expected, xx)
}

test_minor_allele_depth_dropall <- function() {
    scanBamParam <- .min_minor_allele_depth()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=FALSE,
                              min_nucleotide_depth=0L,
                              min_minor_allele_depth=2L)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam))
    expected <- .tinysam_empty_df(scanBamParam)[c("seqnames", "pos", "count", "which_label")]
    ##print(xx);str(xx);print(expected);str(expected)
    checkIdentical(expected, xx)
}

## MIN NUC DEPTH

test_min_nuc_depth_polyread <- function() {
    ## should be implemented in a way that this won't fail while
    ## the other min_nuc_depth tests pass, but since code is still
    ## highly volatile, make sure code is robust to those changes
    scanBamParam <- .min_nuc_depth_poly()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=FALSE,
                              min_nucleotide_depth=2L)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos=  1L,
                     count=2L,
                      which_label=.mwls(scanBamParam, length(seqnames)))
    ##print(xx);print(expected)
    checkIdentical(expected, xx)
}

test_min_nuc_depth_dropnone <- function() {
    ## min_nuc_depth is one alignment with single 'A'
    scanBamParam <- .min_nuc_depth_mono()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=FALSE,
                              min_nucleotide_depth=0L)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos=  1L,
                     count=1L,
                      which_label=.mwls(scanBamParam, length(seqnames)))
    ##print(xx);print(expected)
    checkIdentical(expected, xx)
}
test_min_nuc_depth_dropall <- function() {
    ## min_nuc_depth is one alignment with single 'A'
    scanBamParam <- .min_nuc_depth_mono()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=FALSE,
                              min_nucleotide_depth=2L)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam))
    expected <- .tinysam_empty_df(scanBamParam)[c("seqnames", "pos", "count", "which_label")]
    ##print(xx);str(xx);print(expected);str(expected)
    checkIdentical(expected, xx)
}

## INCLUDE DELETIONS

test_include_dels_dist_nucs <- function() {
    scanBamParam <- .include_deletions()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=TRUE,
                              ignore_query_Ns=FALSE,
                              include_deletions=TRUE)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos=       c(1L, 2L, 3L),
                     nucleotide=c("A", "-", "A"),
                     count=     c(1L, 1L, 1L),
                      which_label=.mwls(scanBamParam, length(seqnames)))
    ##print(xx);print(expected)
    checkIdentical(expected, xx)
}
test_include_dels_no_dist_nucs <- function() {
    scanBamParam <- .include_deletions()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=FALSE,
                              ignore_query_Ns=FALSE,
                              include_deletions=TRUE)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos=  c(1L, 2L, 3L),
                     count=c(1L, 1L, 1L),
                      which_label=.mwls(scanBamParam, length(seqnames)))
    ##print(xx);print(expected)
    checkIdentical(expected, xx)
}
test_no_include_dels_no_dist_nucs <- function() {
    scanBamParam <- .include_deletions()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=FALSE,
                              ignore_query_Ns=FALSE,
                              include_deletions=FALSE)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos=  c(1L, 3L),
                     count=c(1L, 1L),
                      which_label=.mwls(scanBamParam, length(seqnames)))
    ##print(xx);print(expected)
    checkIdentical(expected, xx)
}

## REF SKIP

test_ref_skip <- function() {
    ## sanity check that an N in a cigar means don't tally it--ever!
    ## (regardless of what include_deletions says)
    scanBamParam <- .ref_skip()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=FALSE,
                              ignore_query_Ns=FALSE,
                              include_deletions=TRUE)
    include_dels <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=FALSE,
                              ignore_query_Ns=FALSE,
                              include_deletions=FALSE)
    exclude_dels <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    ##print(include_dels);print(exclude_dels)
    seqnames <- space(bamWhich(scanBamParam))
    expected <- .tadf(seqnames=seqnames,
                     pos=  c(1L, 3L),
                     count=c(1L, 1L),
                      which_label=.mwls(scanBamParam, length(seqnames)))
    ##print(expected)
    checkIdentical(expected, include_dels)
    checkIdentical(include_dels, exclude_dels)
}

## SOFT- AND HARD-CLIPPING

test_ext_cigar_clipping <- function() {
    ## keep in as a sanity check that insert/filter code doesn't mess clips
    scanBamParam <- .ext_cigar_hard()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=FALSE,
                              ignore_query_Ns=FALSE,
                              include_deletions=TRUE)
    hard_res <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    ##print("hard clipping tACttGTN")
    ##print(hard_res)
    scanBamParam <- .ext_cigar_soft()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=FALSE,
                              ignore_query_Ns=FALSE,
                              include_deletions=TRUE)
    soft_res <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    ##print("soft clipping ACGTN")
    ##print(soft_res)
    soft_res$seqnames <- NULL
    soft_res$which_label <- NULL
    hard_res$seqnames <- NULL
    hard_res$which_label <- NULL

    checkIdentical(hard_res, soft_res)
}

## MAX DEPTH

test_max_depth_250 <- function() {
    ## max_depth has 3 alignments
    ## Important case because weird edge cases in samtools (as indicated in
    ## PileupBuffer code)
    scanBamParam <- .max_depth()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=FALSE,
                              ignore_query_Ns=FALSE,
                              max_depth=250L)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos=  1L,
                     count=3L,
                      which_label=.mwls(scanBamParam, length(seqnames)))
    checkIdentical(expected, xx)
}
test_max_depth_2 <- function() {
    ## max_depth has 3 alignments
    ## Important case because weird edge cases in samtools (as indicated in
    ## PileupBuffer code)
    scanBamParam <- .max_depth()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=FALSE,
                              ignore_query_Ns=FALSE,
                              max_depth=2L)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos=  1L,
                     count=2L,
                      which_label=.mwls(scanBamParam, length(seqnames)))
    checkIdentical(expected, xx)
}
test_max_depth_1 <- function() {
    ## max_depth has 3 alignments
    ## Important case because weird edge cases in samtools (as indicated in
    ## PileupBuffer code)
    scanBamParam <- .max_depth()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=FALSE,
                              ignore_query_Ns=FALSE,
                              max_depth=1L)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos=  1L,
                     count=1L,
                      which_label=.mwls(scanBamParam, length(seqnames)))
    checkIdentical(expected, xx)
}
test_max_depth_0_exception <- function() {
    ## max_depth has 3 alignments
    scanBamParam <- .max_depth()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=FALSE,
                              ignore_query_Ns=FALSE,
                              max_depth=0L)
    ##obs <- tryCatch(pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    checkException(f <- tryCatch(pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)), silent=TRUE)
}

## MIN MAPQ

test_min_mapq_dropnone <- function() {
    scanBamParam <- .min_mapq()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=FALSE,
                              ignore_query_Ns=FALSE,
                              min_mapq=5L)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos=  1L,
                     count=3L,
                      which_label=.mwls(scanBamParam, length(seqnames)))
    ##print(xx);print(expected)
    checkIdentical(expected, xx)
}
test_min_mapq_dropsome <- function() {
    scanBamParam <- .min_mapq()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=FALSE,
                              ignore_query_Ns=FALSE,
                              min_mapq=15L)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos=  1L,
                     count=2L,
                      which_label=.mwls(scanBamParam, length(seqnames)))
    ##print(xx);print(expected)
    checkIdentical(expected, xx)
}
test_min_mapq_dropall <- function() {
    scanBamParam <- .min_mapq()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=FALSE,
                              ignore_query_Ns=FALSE,
                              min_mapq=31L)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam))
    expected <- .tinysam_empty_df(scanBamParam)[c("seqnames", "pos", "count", "which_label")]
    ##print(xx);str(xx);print(expected);str(expected)
    checkIdentical(expected, xx)
}

## MIN BASEQUAL

test_min_basequal_dropnone <- function() {
    ## min_base_qual QUAL values correspond to 0, 10, 20
    scanBamParam <- .min_base_qual()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=FALSE,
                              ignore_query_Ns=FALSE,
                              min_base_quality=0L)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos=  c(1L, 2L, 3L),
                     count=c(1L, 1L, 1L),
                      which_label=.mwls(scanBamParam, length(seqnames)))
    ##print(xx);print(expected)
    checkIdentical(expected, xx)
}
test_min_basequal_dropsome <- function() {
    ## min_base_qual QUAL values correspond to 0, 10, 20
    scanBamParam <- .min_base_qual()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=FALSE,
                              ignore_query_Ns=FALSE,
                              min_base_quality=20L)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos=  3L,
                     count=1L,
                      which_label=.mwls(scanBamParam, length(seqnames)))
    ##print(xx);print(expected)
    checkIdentical(expected, xx)
}
test_min_basequal_dropall <- function() {
    ## min_base_qual QUAL values correspond to 0, 10, 20
    scanBamParam <- .min_base_qual()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=FALSE,
                              ignore_query_Ns=FALSE,
                              min_base_quality=21L)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam))
    expected <- .tinysam_empty_df(scanBamParam)[c("seqnames", "pos", "count", "which_label")]
    
    ##print(xx);str(xx);print(expected);str(expected)
    checkIdentical(expected, xx)
}

## DISTINGUISH STRANDS

test_distinguish_strands_plus <- function() {
    scanBamParam <- .seq1()
    pileupParam <- PileupParam(
                              distinguish_strands=TRUE,
                              distinguish_nucleotides=FALSE,
                              ignore_query_Ns=TRUE)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos=       seq(1L, 5L, 1L),
                     strand=   rep("+", 5),
                     count=     rep(1L, 5L),
                      which_label=.mwls(scanBamParam, length(seqnames)))
    checkIdentical(expected, xx)
}
test_distinguish_strands_minus <- function() {
    scanBamParam <- .seq2()
    pileupParam <- PileupParam(
                              distinguish_strands=TRUE,
                              distinguish_nucleotides=FALSE,
                              ignore_query_Ns=FALSE)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos=       seq(1L, 5L, 1L),
                     strand=   rep("-", 5),
                     count=     rep(1L, 5L),
                      which_label=.mwls(scanBamParam, length(seqnames)))
    checkIdentical(expected, xx)
}

## DISTINGUISH NUCLEOTIDES

test_distinguish_nucleotides_include_Ns <- function() {
    scanBamParam <- .seq2()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=TRUE,
                              ignore_query_Ns=FALSE)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos=       seq(1L, 5L, 1L),
                     nucleotide=c("A", "C", "N", "T", "A"),
                     count=     rep(1L, 5L),
                      which_label=.mwls(scanBamParam, length(seqnames)))
    checkIdentical(expected, xx)
}

test_distinguish_nucleotides_ignore_Ns <- function() {
    scanBamParam <- .seq2()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=TRUE,
                              ignore_query_Ns=TRUE)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos=       c(1L, 2L, 4L, 5L),
                     nucleotide=c("A", "C", "T", "A"),
                     count=     rep(1L, 4L),
                      which_label=.mwls(scanBamParam, length(seqnames)))
    checkIdentical(expected, xx)
}

test_distinguish_nothing <- function() {
    scanBamParam <- .seq1()
    pileupParam <- PileupParam(
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=FALSE,
                              ignore_query_Ns=FALSE)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos=       seq(1L, 5L, 1L),
                     count=     rep(1L, 5L),
                      which_label=.mwls(scanBamParam, length(seqnames)))
    checkIdentical(expected, xx)
}

## SCHEMA

test_schema_default_params <- function() { # distinguish all
    scanBamParam <- .seq1()
    pileupParam <- PileupParam()
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos = rep(0L, 5),
                     strand = rep("*", 5),
                     nucleotide = rep("A", 5),
                     count = rep(0L, 5),
                      which_label=.mwls(scanBamParam, length(seqnames)))
    checkIdentical(sapply(expected, names), sapply(xx, names))
    checkIdentical(sapply(expected, length), sapply(xx, length))
    checkIdentical(sapply(expected, class), sapply(xx, class))
}

test_schema_distinguish_single_field <- function() {
    scanBamParam <- .seq1()
    pileupParam <- PileupParam( distinguish_nucleotides=FALSE,
                              distinguish_strands=TRUE)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos = rep(0L, 5),
                     strand = rep("*", 5),
                     count = rep(0L, 5),
                      which_label=.mwls(scanBamParam, length(seqnames)))
    checkIdentical(sapply(expected, names), sapply(xx, names))
    checkIdentical(sapply(expected, length), sapply(xx, length))
    checkIdentical(sapply(expected, class), sapply(xx, class))
}

test_schema_distinguish_nothing <- function() {
    scanBamParam <- .seq1()
    pileupParam <- PileupParam( distinguish_nucleotides=FALSE,
                              distinguish_strands=FALSE)
    xx <- pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
    seqnames <- space(bamWhich(scanBamParam)) 
    expected <- .tadf(seqnames=seqnames,
                     pos = rep(0L, 5),
                     count = rep(0L, 5),
                      which_label=.mwls(scanBamParam, length(seqnames)))
    checkIdentical(sapply(expected, names), sapply(xx, names))
    checkIdentical(sapply(expected, length), sapply(xx, length))
    checkIdentical(sapply(expected, class), sapply(xx, class))
}

test_schemabuilder_default_params <- function() {
    pileupParam <- PileupParam()
    schema <- Rsamtools:::.schemaBuilder(pileupParam)
    expected <- list(c("strand", "nucleotide"),
                     list(c("+", "-"), c("A", "C", "G", "T", "-")))
    checkIdentical(expected, schema)
}

test_schemabuilder_distinguish_nothing <- function() {
    pileupParam <- PileupParam(distinguish_nucleotides=FALSE, # non-default
                              distinguish_strands=FALSE) # non-default
    schema <- Rsamtools:::.schemaBuilder(pileupParam)
    expected <- list(c("strand", "nucleotide"),
                     list("", ""))
    checkIdentical(expected, schema)
}

## REVERSE COMPLEMENT

test_reverseComplement_warning <- function() {
    scanBamParam <- ScanBamParam(reverseComplement=TRUE)
    pileupParam <- PileupParam()
    obs <- tryCatch(pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam), warning=conditionMessage)
    expectedWarning <-
        "'reverseComplement' parameter in pileup ScanBamParam will be ignored"
    checkIdentical(expectedWarning, obs)
}
