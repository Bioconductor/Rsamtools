test_readBAMasGappedAlignments <- function()
{
    fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
    which <- RangesList(seq1=IRanges(1, width=100))
    result <- readBAMasGappedAlignments(fl, which=which)
    checkTrue(validObject(result))
}

test_readBAMasGappedAlignments_length0 <- function()
{
    fl <- system.file("extdata", "ex1.bam", package="Rsamtools")

    which <- RangesList(seq1=IRanges(100000, width=100))
    result <- readBAMasGappedAlignments(fl, which=which)
    checkTrue(validObject(result))

    which <- RangesList(seq1=IRanges(c(1, 100000), width=100))
    result <- readBAMasGappedAlignments(fl, which=which)
    checkTrue(validObject(result))
}
