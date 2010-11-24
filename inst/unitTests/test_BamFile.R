library(Rsamtools); library(RUnit)
fl <- system.file("extdata", "ex1.bam", package="Rsamtools")

test_BamFile_openclose <- function()
{
    .normalizePath <- Rsamtools:::.normalizePath
    bf <- open(BamFile(fl))
    checkTrue(isOpen(bf))
    checkIdentical(.normalizePath(fl), path(bf))
    checkIdentical(.normalizePath(fl), index(bf))
    close(bf)
    checkTrue(!isOpen(bf))
    checkException(close(bf), silent=TRUE)
    bf <- open(bf)                   # open a closed BamFile
    checkTrue(isOpen(bf))
    bf1 <- open(bf)                  # open (clone?) an open BamFile
    checkTrue(isOpen(bf1))
    checkTrue(identical(bf$.extptr, bf1$.extptr))
}
