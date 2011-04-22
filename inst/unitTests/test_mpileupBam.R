if (FALSE) {
    library(Rsamtools); library(RUnit)
    fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
    MpileupParam <- Rsamtools:::.MpileupParam
    mpileupBam <- Rsamtools:::.mpileupBam

    fls <- list(open(BamFile(fl)), open(BamFile(fl)))

    callback <-
        function(x)
        {
            ## x is a list(base_calls, quality)
            ## base_calls is a base x cycle x file array
            ## quality (will be) a quality x cycle x file array
            ## e.g., reduce to per-cycle information content
            apply(x[[1]], 3, function(y) {
                y <- y[c(1, 2, 4, 8) + 1,]      # A, C, G, T
                y <- y + 1L                     # continuity
                cvg <- colSums(y)
                p <- y / cvg[col(y)]
                h <- -colSums(p * log(p))
                ifelse(cvg == 4L, NA, h)
            })
        }

    param <- MpileupParam(which=GRanges("seq1", IRanges(1000,2000)))
    res <- mpileupBam(fls, callback, param=param)

    param <- MpileupParam(flag=scanBamFlag(isProperPair=TRUE),
                          minDepth=20L,
                          which=GRanges("seq1", IRanges(1000, 2000)))
    res <- mpileupBam(fls, callback, param=param)
}
