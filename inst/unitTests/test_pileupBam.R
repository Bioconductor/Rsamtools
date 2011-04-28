library(Rsamtools); library(RUnit)
fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
PileupParam <- Rsamtools:::.PileupParam
pileupBam <- Rsamtools:::.pileupBam

if (interactive()) {
    fls <- list(open(BamFile(fl)), open(BamFile(fl)))

    callback <-
        function(x)
        {
            ## x is a list(pos, base_calls, quality)
            ## pos is integer() of length equal to positions passing filter
            ## base_calls is a base x file x cycle array
            ## quality (will be) a quality x cycle x file array
            ## e.g., reduce to per-cycle information content
            list(x[[1]], apply(x[[2]], 2, function(y) {
                y <- y[c(1, 2, 4, 8) + 1,]      # A, C, G, T
                y <- y + 1L                     # continuity
                cvg <- colSums(y)
                p <- y / cvg[col(y)]
                h <- -colSums(p * log(p))
                ifelse(cvg == 4L, NA, h)
            }))
        }
    param <- PileupParam(which=GRanges("seq1", IRanges(1000,2000)))
    res <- pileupBam(fls, callback, param=param)

    param <- PileupParam(flag=scanBamFlag(isProperPair=TRUE),
                          minDepth=20L,
                          which=GRanges("seq1", IRanges(1000, 2000)))
    res <- pileupBam(fls, callback, param=param)

    ## a bigger example
    fls0 <-
        list.files("/home/mtmorgan/a/bioC/Courses/Seattle-Dec-2010/content/SeattleIntro2010/NagalakshmiEtAl/aln",
                   pattern=".*sorted.bam$", full=TRUE)

    fls <- lapply(fls0, function(fl) open(BamFile(fl)))
    scanBamHeader(fls[[1]])[["targets"]]
    irng <- successiveIRanges(rep(200000, 5), from=1)
    param <- PileupParam(which=GRanges("chrIV", irng),
                         minDepth=40L)
    res <- pileupBam(fls, function(x) length(x[[1]]), param=param)
}

