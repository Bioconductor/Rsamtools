## library(Rsamtools); library(RUnit)
fl <- system.file("extdata", "ex1.bam", package="Rsamtools")

test_applyPileups_byRange <- function()
{
    fls <- PileupFiles(c(fl, fl))
    fun <- function(x) x[["seqnames"]]
    which <- GRanges(c("seq1", "seq2"), IRanges(c(1000, 1000), 2000))
    param <- ApplyPileupsParam(which=which)
    res <- applyPileups(fls, fun, param=param)
    exp <- list(structure(570L, .Names = "seq1"),
                structure(568L, .Names = "seq2"))
    checkIdentical(exp, res)

    fls <- PileupFiles(c(fl, fl), param=param)
    res <- applyPileups(fls, fun)
    checkIdentical(exp, res)
}

test_applyPileups_byRange_yieldAll <- function()
{
    fls <- PileupFiles(c(fl, fl))
    which <- GRanges(c("seq1", "seq2"),
                     IRanges(c(1000, 1000), c(2000, 1500)))
    param <- ApplyPileupsParam(which=which, yieldAll=TRUE)
    res <- applyPileups(fls, function(x) x$seqnames, param=param)
    exp <- list(structure(1001L, .Names = "seq1"),
                structure(501L, .Names = "seq2"))
    checkIdentical(exp, res)
    res <- applyPileups(fls, function(x) x$pos, param=param)
    checkIdentical(list(1000:2000, 1000:1500), res)
}

test_applyPileups_byPosition <- function()
{
    fls <- PileupFiles(c(fl, fl))
    fun <- function(x) x[["seqnames"]]
    which <- GRanges(c("seq1", "seq2"), IRanges(c(1000, 1000), 2000))
    param <-
        ApplyPileupsParam(which=which, yieldSize=2000L, yieldBy="position")
    res <- applyPileups(fls, fun, param=param)
    exp <- list(structure(c(570L, 568L), .Names = c("seq1", "seq2")))
    checkIdentical(exp, res)

    res <- applyPileups(fls, function(x) sum(x$seq[,1,]), param=param)
    checkIdentical(list(41165L), res)

    param <-
        ApplyPileupsParam(which=which, yieldSize=500L, yieldBy="position")
    res <- applyPileups(fls, fun, param=param)
    exp <- list(structure(500L, .Names = "seq1"),
                structure(c(70L, 430L), .Names = c("seq1", "seq2")),
                structure(138L, .Names = "seq2"))
    checkIdentical(exp, res)

    res <- applyPileups(fls, function(x) sum(x$seq[,1,]), param=param)
    checkIdentical(list(18132L, 20125L, 2908L), res)
}

test_applyPileups_byPosition_yieldAll <- function() {
    fls <- PileupFiles(c(fl, fl))
    fun <- function(x) colSums(x$seq[,1,,drop=FALSE])
    which <- GRanges("seq1", IRanges(1000, 1100))

    param <- ApplyPileupsParam(which=which, yieldSize=10L,
                               yieldBy="position", yieldAll=TRUE)
    res0 <- applyPileups(fls, fun, param=param)
    checkIdentical(11L, length(res0))

    param <- ApplyPileupsParam(which=which, yieldSize=50L,
                               yieldBy="position", yieldAll=TRUE)
    res <- applyPileups(fls, fun, param=param)
    checkIdentical(3L, length(res))
    checkIdentical(unlist(res0), unlist(res))

    param <- ApplyPileupsParam(which=which, yieldSize=100L,
                               yieldBy="position", yieldAll=TRUE)
    res <- applyPileups(fls, fun, param=param)
    checkIdentical(2L, length(res))
    checkIdentical(unlist(res0), unlist(res))

    param <- ApplyPileupsParam(which=which, yieldSize=101L,
                               yieldBy="position", yieldAll=TRUE)
    res <- applyPileups(fls, fun, param=param)
    checkIdentical(unlist(res0), unlist(res))
    checkIdentical(1L, length(res))

    param <- ApplyPileupsParam(which=which, yieldSize=200L,
                               yieldBy="position", yieldAll=TRUE)
    res <- applyPileups(fls, fun, param=param)
    checkIdentical(1L, length(res))
    checkIdentical(unlist(res0), unlist(res))

    param <- ApplyPileupsParam(which=GRanges("seq1",
                               IRanges(seq(10, 110, by=10), width=1), "+"),
                               yieldBy="position", yieldAll=TRUE, yieldSize=5)
    seqnames <- applyPileups(fls, function(x) x$seqnames, param=param)
    exp <- list(structure(5L, .Names = "seq1"),
                structure(5L, .Names = "seq1"),
                structure(1L, .Names = "seq1"))
    checkIdentical(exp, seqnames)
    pos <- applyPileups(fls, function(x) x$pos, param=param)
    checkIdentical(c(5L, 5L, 1L), sapply(pos, length))
    checkIdentical(seq(10L, 110L, 10L), unlist(pos))
    fun <- function(x) as.vector(colSums(x$seq[,1,,drop=FALSE]))
    seq <- applyPileups(fls, fun, param=param)
    checkIdentical(list(c(4, 8, 13, 12, 10), c(13, 15, 17, 9, 3), 5), seq)
}

test_applyPileups_what <- function() {
    fls <- PileupFiles(c(fl, fl))
    which <- GRanges("seq1", IRanges(1000, 1999))

    param <- ApplyPileupsParam(which=which, yieldAll=TRUE)
    obs <- applyPileups(fls, function(x) sapply(x, length), param=param)[[1]]
    exp <- c(seqnames=1L, pos=1000L, seq=2L * 5L * 1000L,
             qual=2L * 94L * 1000L)
    checkIdentical(obs, exp)

    param <- ApplyPileupsParam(which=which, yieldAll=TRUE, what=character())
    obs <- applyPileups(fls, function(x) x$pos, param=param)[[1]]
    checkIdentical(999L + 1:1000, obs)

    param <- ApplyPileupsParam(which=which, yieldAll=TRUE, what="seq")
    obs <- applyPileups(fls, function(x) sapply(x, length), param=param)[[1]]
    checkIdentical(obs, exp[1:3])

    param <- ApplyPileupsParam(which=which, yieldAll=TRUE, what="qual")
    obs <- applyPileups(fls, function(x) sapply(x, length), param=param)[[1]]
    checkIdentical(obs, exp[c(1:2, 4)])
}

test_applyPileups_memoryleak_warning <- function() {
    ## failing to complete an iterator causes a warning; this is
    ## corrected in C code
    opts <- options(warn=2)
    on.exit(options(opts))
    fls <- PileupFiles(c(fl, fl))
    param <- ApplyPileupsParam(which=GRanges("seq1", IRanges(1000, 1499)),
                               yieldAll=TRUE)
    obs <- applyPileups(fls, function(x) NULL, param=param)
    checkIdentical(list(NULL), obs)
}

test_applyPileups_NULL_space <- function() {
    checkException(applyPileups(PileupFiles(fl), identity),
                   silent=TRUE)
}

test_applyPileups_skip_inserts <- function() {
    fl <- system.file("unitTests", "cases", "plp_refskip.bam", 
                      package="Rsamtools")
    param <- ApplyPileupsParam(which=GRanges("seq1", IRanges(1, 1000000)),
                               what="seq")
    res <-
        applyPileups(PileupFiles(fl), function(x) rowSums(x$seq[,1,]),
                     param=param)

    obs <- structure(c(19, 12, 18, 26, 0),
                     .Names = c("A", "C", "G", "T", "N"))
    checkIdentical(obs, res[[1]])

    res <-
        applyPileups(PileupFiles(fl), function(x) length(x$pos),
                     param=param)
    checkIdentical(75L, res[[1]])

    ## yieldAll=TRUE
    param <- ApplyPileupsParam(which=GRanges("seq1", IRanges(1, 1000000)),
                               yieldAll=TRUE, what="seq")
    obs <- structure(c(19, 12, 18, 26, 0),
                     .Names = c("A", "C", "G", "T", "N"))
    res <-
        applyPileups(PileupFiles(fl), function(x) length(x$pos),
                     param=param)
    checkIdentical(width(plpWhich(param)), res[[1]])
}
