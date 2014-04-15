fl <- system.file("extdata", "ex1.bam", package="Rsamtools")

test_BamSampler <- function()
{
    bf <- BamSampler(fl, yieldSize=1000)
    set.seed(123L); obs1 <- scanBam(bf)
    checkIdentical(yieldSize(bf), length(obs1[[1]][[1]]))
    checkIdentical(741747L, sum(obs1[[1]][["pos"]], na.rm=TRUE))

    set.seed(123L); obs2 <- scanBam(bf)
    checkIdentical(obs1, obs2)

    obs2 <- scanBam(bf)
    checkTrue(!identical(obs1, obs2))

    bf <- BamSampler(fl, yieldSize=10000)
    obs <- scanBam(bf)
    exp <- list(BamSamplerStatistics=c(yieldSize=10000L,
                  totalRead=3307L, yield=3307L))
    checkIdentical(exp, attributes(obs))
    attributes(obs) <- NULL
    checkIdentical(scanBam(fl), obs)
}

test_BamSampler_param <- function()
{
    bf <- BamSampler(fl, yieldSize=1000)
    param <- ScanBamParam(what="strand", flag=scanBamFlag(
                                           isUnmappedQuery=FALSE,
                                           isMinusStrand=FALSE))
    obs <- scanBam(bf, param=param)[[1]][["strand"]]
    checkTrue(all(obs == "+"))
    checkIdentical(yieldSize(bf), length(obs))

    bf <- BamSampler(fl, yieldSize=10000)
    obs <- scanBam(bf, param=param)[[1]][["strand"]]
    checkIdentical(1647L, length(obs))
}
