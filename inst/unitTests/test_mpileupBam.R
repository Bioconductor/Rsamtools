library(Rsamtools); library(RUnit)
fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
mpileupBam <- Rsamtools:::.mpileupBam

test_mpileupBam <-function()
{
    param <- MpileupParam(which=GRanges("seq1", IRanges(1000,2000)))
    fls <- list(open(BamFile(fl)), open(BamFile(fl)))
    mpileupBam(fls, param=param)
}
