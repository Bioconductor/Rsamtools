library(Rsamtools); library(RUnit)
fl <- system.file("extdata", "ex1.bam", package="Rsamtools")

test_BamFile_openclose <- function()
{
    bf <- openBam(fl)
    checkTrue(isOpen(bf))
    checkIdentical(fl, bamPath(bf))
    checkIdentical(fl, bamIndex(bf))
    res <- close(bf)
    checkTrue(!isOpen(bf))
    checkException(close(bf), silent=TRUE)
    bf <- openBam(bf)                   # open a closed BamFile
    checkTrue(isOpen(bf))
    bf1 <- openBam(bf)                  # open (clone?) an open BamFile
    checkTrue(isOpen(bf1))
    checkTrue(!identical(bf$.extptr, bf1$.extptr))
}
