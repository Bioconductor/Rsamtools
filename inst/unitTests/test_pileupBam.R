## library(Rsamtools); library(RUnit)
fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
PileupParam <- Rsamtools:::.PileupParam
pileupBam <- Rsamtools:::.pileupBam

test_pileupBam_byRange <- function()
{
    fls <- list(open(BamFile(fl)), open(BamFile(fl)))
    fun <- function(x) x[["seqnames"]]
    which <- GRanges(c("seq1", "seq2"), IRanges(c(1000, 1000), 2000))
    param <- PileupParam(which=which)
    res <- pileupBam(fls, fun, param=param)
    exp <- list(structure(569L, .Names = "seq1"),
                structure(567L, .Names = "seq2"))
    checkIdentical(exp, res)
}

test_pileupBam_byPosition <- function()
{
    fls <- list(open(BamFile(fl)), open(BamFile(fl)))
    fun <- function(x) x[["seqnames"]]
    which <- GRanges(c("seq1", "seq2"), IRanges(c(1000, 1000), 2000))
    param <-
        PileupParam(which=which, yieldSize=2000L, yieldBy="position")
    res <- pileupBam(fls, fun, param=param)
    exp <- list(structure(c(569L, 567L), .Names = c("seq1", "seq2")))
    checkIdentical(exp, res)

    param <-
        PileupParam(which=which, yieldSize=500L, yieldBy="position")
    res <- pileupBam(fls, fun, param=param)
    exp <- list(structure(500L, .Names = "seq1"),
                structure(c(69L, 431L), .Names = c("seq1", "seq2")),
                structure(136L, .Names = "seq2"))
    checkIdentical(exp, res)
}
