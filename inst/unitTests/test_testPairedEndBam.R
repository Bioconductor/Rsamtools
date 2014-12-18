fl <- system.file("extdata", "ex1.bam", package="Rsamtools")

test_testPairedEndBam <- function() {
    checkTrue(testPairedEndBam(BamFile(fl)))
    checkTrue(testPairedEndBam(fl))
}
