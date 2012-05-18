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

test_BamFile_yield <- function()
{
    bf <- open(BamFile(fl, yieldSize=1000))
    it <- integer()
    while(length(res <- scanBam(bf)[[1]][[1]]))
        it <- append(it, length(res))
    close(bf)
    checkIdentical(c(1000L, 1000L, 1000L, 307L), it)

    open(bf)
    rng <- GRanges(c("seq1", "seq2"), IRanges(1, c(1575, 1584)))
    param <- ScanBamParam(which=rng)
    checkException(scanBam(bf, param=param), silent=TRUE)
    close(bf)
}

test_BamFileList_yield <- function()
{
    bfl <- unname(BamFileList(c(fl, fl), yieldSize=100))
    checkIdentical(c(100L, 100L), sapply(bfl, yieldSize))

    bfl <- unname(BamFileList(c(fl, fl)))
    checkIdentical(c(NA_integer_, NA_integer_), sapply(bfl, yieldSize))

    ## yieldSize, even implicit on BamFile wins out
    bfl <- BamFileList(BamFile(fl, yieldSize=10), BamFile(fl, yieldSize=20),
                       yieldSize=100)
    checkIdentical(c(10L, 20L), sapply(bfl, yieldSize))

    bfl <- BamFileList(BamFile(fl), BamFile(fl), yieldSize=100)
    checkIdentical(c(NA_integer_, NA_integer_), sapply(bfl, yieldSize))
}
