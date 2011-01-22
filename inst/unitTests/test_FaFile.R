fl <- system.file("extdata", "ce2dict1.fa", package="Rsamtools")

test_FaFile_openclose <- function()
{
    fa <- FaFile(fl)
    checkIdentical(FALSE, isOpen(fa))
    checkIdentical(TRUE, isOpen(open(fa)))
    checkIdentical(FALSE, isOpen(close(fa)))
}

test_FaFile_scanFaIndex <- function()
{
    .checkIdx <- function(idx) {
        checkTrue(is(idx, "GRanges"))
        checkIdentical(5L, length(idx))
        checkIdentical(116L, sum(width(idx)))
    }
    fa <- FaFile(fl)
    .checkIdx(scanFaIndex(fa))
    .checkIdx(scanFaIndex(fl))
}

test_FaFile_count <- function()
{
    checkIdentical(5L, countFa(open(FaFile(fl))))
    checkIdentical(5L, countFa(fl))
}

test_FaFile_scanFa <- function()
{
    .checkRes <- function(res) {
        checkTrue(is(res, "DNAStringSet"))
        checkIdentical(5L, length(idx))
        checkIdentical(116L, sum(width(idx)))
    }
    fa <- open(FaFile(fl))
    idx <- scanFaIndex(fa)[1:5]
    .checkRes(scanFa(fa, idx))
    .checkRes(scanFa(fl, idx))

    ## scanFa,*,missing-methods
    checkTrue(validObject(scanFa(fa)))
    checkTrue(validObject(scanFa(fl)))
}
