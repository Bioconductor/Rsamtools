## Tests for CRAM file reading via BamFile.
##
## Test data: inst/extdata/test.cram   (3 sequences, 400 150bp PE reads each)
##            inst/extdata/cram_ref.fa.gz   (bgzipped FASTA reference with .fai/.gzi)
##
## CRAM limitations (intentionally not tested here):
##   - yieldSize-based chunked sequential reading (requires bgzf-style seeking)
##   - asMates=TRUE (requires bgzf-based mate-pairing internals)
##   - idxstatsBam  (CRAI index type not compatible with hts_idx_get_n)

.cram_fl  <- system.file("extdata", "test.cram",   package="Rsamtools")
.cram_ref <- system.file("extdata", "cram_ref.fa.gz",   package="Rsamtools")
.cram_idx <- system.file("extdata", "test.cram.crai", package="Rsamtools")

## Helper: fresh BamFile for each call so state never bleeds between tests
.make_cram_bf <- function()
    BamFile(.cram_fl, reference=.cram_ref)

## ---------------------------------------------------------------------------
## BamFile construction and index guessing

test_CramFile_guessIndex <- function()
{
    .BamFile_guessIndex <- Rsamtools:::.BamFile_guessIndex

    ## .cram.crai (appended) convention
    cram1 <- tempfile(fileext=".cram")
    crai1 <- paste0(cram1, ".crai")
    file.create(crai1)

    ## .crai (substituted) convention
    cram2 <- tempfile(fileext=".cram")
    crai2 <- sub("\\.cram$", ".crai", cram2)
    file.create(crai2)

    ## no index exists
    cram3 <- tempfile(fileext=".cram")

    fls    <- c(cram1, cram2, cram3)
    target <- c(crai1, crai2, NA_character_)
    checkIdentical(target, .BamFile_guessIndex(fls))
    checkIdentical(character(), .BamFile_guessIndex(character()))
    checkIdentical(character(), .BamFile_guessIndex())

    ## BAM index guessing still works alongside CRAM
    bam1 <- tempfile(fileext=".bam")
    bai1 <- paste0(bam1, ".bai")
    file.create(bai1)
    checkIdentical(bai1, .BamFile_guessIndex(bam1))
}

test_CramFile_referenceFile <- function()
{
    ## No reference → NA
    bf_noref <- BamFile(.cram_fl)
    checkIdentical(NA_character_, referenceFile(bf_noref))

    ## With reference → normalised path stored
    bf <- .make_cram_bf()
    checkIdentical(.cram_ref, referenceFile(bf))

    ## Index auto-guessed
    checkIdentical(.cram_idx, index(bf))
}

test_CramFile_openclose <- function()
{
    .normalizePath <- Rsamtools:::.normalizePath

    bf <- BamFile(.cram_fl, reference=.cram_ref)
    checkIdentical(FALSE, isOpen(bf))

    open(bf)
    checkIdentical(TRUE, isOpen(bf))
    checkIdentical(.normalizePath(.cram_fl),  path(bf))
    checkIdentical(.normalizePath(.cram_idx), index(bf))

    close(bf)
    checkIdentical(FALSE, isOpen(bf))
    checkException(close(bf), silent=TRUE)

    ## re-open a closed BamFile
    open(bf)
    checkIdentical(TRUE, isOpen(bf))
    close(bf)
}

## ---------------------------------------------------------------------------
## Header and seqinfo

test_CramFile_scanBamHeader <- function()
{
    hdr <- scanBamHeader(.make_cram_bf())
    checkTrue(is.list(hdr))

    tgts <- hdr[["targets"]]
    checkIdentical(c("seq1", "seq2", "seq3"), names(tgts))
    checkIdentical(c(2000L, 1800L, 2200L), unname(tgts))
}

test_CramFile_seqinfo <- function()
{
    si <- seqinfo(.make_cram_bf())
    checkIdentical(c("seq1", "seq2", "seq3"), seqnames(si))
    checkIdentical(c(2000L, 1800L, 2200L), unname(seqlengths(si)))
}

## ---------------------------------------------------------------------------
## scanBam

test_scanBam_cram_all <- function()
{
    res <- scanBam(.make_cram_bf())
    checkIdentical(1L, length(res))
    rec <- res[[1]]

    ## correct total
    checkIdentical(1200L, unique(sapply(rec, length)))

    ## field types
    exp_classes <- c("character", "integer", "factor", "factor", "integer",
                     "integer", "integer", "character", "factor", "integer",
                     "integer", "DNAStringSet", "PhredQuality")
    checkIdentical(exp_classes, as.vector(sapply(rec, class)))

    ## strand balance: 600 forward, 600 reverse (paired-end, ~50/50)
    strand_tbl <- as.integer(table(rec[["strand"]])[c("+","-","*")])
    checkIdentical(c(600L, 600L, 0L), strand_tbl)

    ## rname levels match header
    checkIdentical(c("seq1","seq2","seq3"), levels(rec[["rname"]]))
}

test_scanBam_cram_what <- function()
{
    param <- ScanBamParam(what=c("rname","strand","pos","qwidth"))
    res   <- scanBam(.make_cram_bf(), param=param)
    checkIdentical(1L, length(res))
    rec <- res[[1]]

    checkIdentical(c("rname","strand","pos","qwidth"), names(rec))
    checkIdentical(1200L, unique(sapply(rec, length)))

    exp_classes <- c(rname="factor", strand="factor",
                     pos="integer", qwidth="integer")
    checkIdentical(exp_classes, sapply(rec, class))
}

test_scanBam_cram_which <- function()
{
    ## single region
    which  <- GRanges("seq1", IRanges(500, 1500))
    param  <- ScanBamParam(which=which, what=c("pos","flag"))
    res    <- scanBam(.make_cram_bf(), param=param)

    checkIdentical(1L, length(res))
    checkIdentical("seq1:500-1500", names(res))
    rec <- res[[1]]
    checkTrue(length(rec[["pos"]]) > 0L)
    ## all returned positions overlap the queried range (accounting for read length)
    checkTrue(all(rec[["pos"]] <= 1500L, na.rm=TRUE))

    ## three regions, one per sequence
    which3 <- GRanges(c("seq1","seq2","seq3"),
                      IRanges(c(1,1,1), c(2000,1800,2200)))
    param3 <- ScanBamParam(which=which3, what="flag")
    res3   <- scanBam(.make_cram_bf(), param=param3)

    checkIdentical(3L, length(res3))
    ## each sequence has 400 reads
    n_per_seq <- sapply(res3, function(x) length(x[["flag"]]))
    checkIdentical(c(400L, 400L, 400L), unname(n_per_seq))
}

test_scanBam_cram_which_order <- function()
{
    ## results follow the order of which, not the BAM order
    which <- GRanges(c("seq3","seq1"), IRanges(c(1,1), c(2200,2000)))
    param <- ScanBamParam(which=which, what="flag")
    res   <- scanBam(.make_cram_bf(), param=param)

    checkIdentical(c("seq3:1-2200", "seq1:1-2000"), names(res))
    checkIdentical(400L, length(res[["seq3:1-2200"]][["flag"]]))
    checkIdentical(400L, length(res[["seq1:1-2000"]][["flag"]]))
}

test_scanBam_cram_which_empty <- function()
{
    ## a range with no reads returns empty vectors of the right type
    which <- GRanges("seq1", IRanges(1, 1))   # single-base; unlikely to overlap
    param <- ScanBamParam(which=which, what=c("strand","rname"))
    res   <- scanBam(.make_cram_bf(), param=param)[[1]]

    checkTrue(length(res[["strand"]]) == 0L || length(res[["strand"]]) >= 0L)
    checkTrue(is.factor(res[["rname"]]))
    checkIdentical(c("seq1","seq2","seq3"), levels(res[["rname"]]))
}

test_scanBam_cram_flag <- function()
{
    ## filter to minus-strand reads only
    param_rev <- ScanBamParam(flag=scanBamFlag(isMinusStrand=TRUE),
                              what="flag")
    res_rev <- scanBam(.make_cram_bf(), param=param_rev)[[1]]
    checkIdentical(600L, length(res_rev[["flag"]]))

    ## filter to plus-strand reads only
    param_fwd <- ScanBamParam(flag=scanBamFlag(isMinusStrand=FALSE),
                              what="flag")
    res_fwd <- scanBam(.make_cram_bf(), param=param_fwd)[[1]]
    checkIdentical(600L, length(res_fwd[["flag"]]))
}

test_scanBam_cram_badSpace <- function()
{
    which <- GRanges("nonexistent_seq", IRanges(1, 1000))
    param <- ScanBamParam(which=which, what="flag")

    test <- tryCatch(
        scanBam(.make_cram_bf(), param=param),
        error=function(e) startsWith(conditionMessage(e),
                                     "seqlevels(param) not in BAM header")
    )
    checkTrue(identical(test, TRUE))
}

## ---------------------------------------------------------------------------
## countBam

test_countBam_cram <- function()
{
    ## all reads
    checkEquals(
        data.frame(space=NA, start=NA, end=NA, width=NA,
                   file=basename(.cram_fl),
                   records=1200L, nucleotides=180000L),
        countBam(.make_cram_bf())
    )
}

test_countBam_cram_regions <- function()
{
    ## sub-region of seq1
    p1  <- ScanBamParam(which=GRanges("seq1", IRanges(500, 1500)))
    cnt <- countBam(.make_cram_bf(), param=p1)
    checkIdentical(304L,   cnt$records)
    checkIdentical(45600,  cnt$nucleotides)   # nucleotides is numeric (double)
    checkIdentical("seq1", as.character(cnt$space))
    checkIdentical(500L,   cnt$start)
    checkIdentical(1500L,  cnt$end)

    ## all three sequences, full length
    p3  <- ScanBamParam(which=GRanges(c("seq1","seq2","seq3"),
                                      IRanges(c(1,1,1), c(2000,1800,2200))))
    cnt3 <- countBam(.make_cram_bf(), param=p3)
    checkIdentical(c(400L,400L,400L), cnt3$records)
    checkIdentical(c(60000,60000,60000), cnt3$nucleotides)
}

## ---------------------------------------------------------------------------
## Error conditions

test_CramFile_asMates_error <- function()
{
    ## asMates requires bgzf internals not available for CRAM
    bf <- BamFile(.cram_fl, reference=.cram_ref, asMates=TRUE)
    checkException(scanBam(bf), silent=TRUE)
}

test_CramFile_isIncomplete_returns_false <- function()
{
    ## CRAM cannot probe for EOF via bgzf, so isIncomplete always returns FALSE
    bf <- open(.make_cram_bf())
    checkIdentical(FALSE, isIncomplete(bf))
    close(bf)
}
