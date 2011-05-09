## library(Rsamtools); library(RUnit)
fl <- system.file("extdata", "ex1.bam", package="Rsamtools")

test_applyPileups_byRange <- function()
{
    fls <- PileupFiles(fl, fl)
    fun <- function(x) x[["seqnames"]]
    which <- GRanges(c("seq1", "seq2"), IRanges(c(1000, 1000), 2000))
    param <- PileupParam(which=which)
    res <- applyPileups(fls, fun, param=param)
    exp <- list(structure(570L, .Names = "seq1"),
                structure(568L, .Names = "seq2"))
    checkIdentical(exp, res)

    fls <- PileupFiles(fl, fl, param=param)
    res <- applyPileups(fls, fun)
    checkIdentical(exp, res)
}

test_applyPileups_byPosition <- function()
{
    fls <- PileupFiles(fl, fl)
    fun <- function(x) x[["seqnames"]]
    which <- GRanges(c("seq1", "seq2"), IRanges(c(1000, 1000), 2000))
    param <-
        PileupParam(which=which, yieldSize=2000L, yieldBy="position")
    res <- applyPileups(fls, fun, param=param)
    exp <- list(structure(c(570L, 568L), .Names = c("seq1", "seq2")))
    checkIdentical(exp, res)

    res <- applyPileups(fls, function(x) sum(x$seq[,1,]), param=param)
    checkIdentical(list(41166L), res)

    param <-
        PileupParam(which=which, yieldSize=500L, yieldBy="position")
    res <- applyPileups(fls, fun, param=param)
    exp <- list(structure(500L, .Names = "seq1"),
                structure(c(70L, 430L), .Names = c("seq1", "seq2")),
                structure(138L, .Names = "seq2"))
    checkIdentical(exp, res)

    res <- applyPileups(fls, function(x) sum(x$seq[,1,]), param=param)
    checkIdentical(list(18132L, 20125L, 2909L), res)
}

test_applyPileups_what <- function() {
    fls <- PileupFiles(fl, fl)
    which <- GRanges("seq1", IRanges(1000, 1999))

    param <- PileupParam(which=which, yieldAll=TRUE)
    obs <- applyPileups(fls, function(x) sapply(x, length), param=param)[[1]]
    exp <- c(seqnames=1L, pos=1000L, seq=32000L, qual=512000L)
    checkIdentical(obs, exp)

    param <- PileupParam(which=which, yieldAll=TRUE, what=character())
    obs <- applyPileups(fls, function(x) x$pos, param=param)[[1]]
    checkIdentical(999L + 1:1000, obs)

    param <- PileupParam(which=which, yieldAll=TRUE, what="seq")
    obs <- applyPileups(fls, function(x) sapply(x, length), param=param)[[1]]
    checkIdentical(obs, exp[1:3])

    param <- PileupParam(which=which, yieldAll=TRUE, what="qual")
    obs <- applyPileups(fls, function(x) sapply(x, length), param=param)[[1]]
    checkIdentical(obs, exp[c(1:2, 4)])
}

test_applyPileups_memoryleak_warning <- function() {
    ## failing to complete an iterator causes a warning; this is
    ## corrected in C code
    opts <- options(warn=2)
    on.exit(options(opts))
    fls <- PileupFiles(fl, fl)
    param <- PileupParam(which=GRanges("seq1", IRanges(1000, 1499)),
                         yieldAll=TRUE)
    obs <- applyPileups(fls, function(x) NULL, param=param)
    checkIdentical(list(NULL), obs)
}

test_applyPileups_NULL_space <- function() {
    checkException(applyPileups(PileupFiles(fl), identity),
                   silent=TRUE)
}
