library(pasillaBamSubset)
readBamSortedPairs <- Rsamtools:::.readBamSortedPairs
fl <- sortBam(untreated3_chr4(), tempfile(), byQname=TRUE)

test_readBamSortedPairs <- function()
{
    bf <- BamFile(fl, index=character(0), yieldSize=100, obeyQname=TRUE)
    galp <- readBamSortedPairs(bf, asGappedAlignmentPairs=TRUE)
    checkTrue(validObject(galp))
    checkTrue(length(galp) == 50)

    gal <- readBamSortedPairs(bf)
    checkTrue(validObject(gal))
    checkTrue(length(gal) == 100)
}

test_readBamSortedPairs_param <- function()
{
    bf <- BamFile(fl, index=character(0), yieldSize=100, obeyQname=TRUE)

    ## valid
    param <- ScanBamParam(tag=("NM"))
    galp <- readBamSortedPairs(bf, param=param, asGappedAlignmentPairs=TRUE)
    checkTrue(names(mcols(left(galp))) == "NM")

    ## empty
    param <- ScanBamParam(tag=("FO"))
    galp <- readBamSortedPairs(bf, param=param, asGappedAlignmentPairs=TRUE)
    checkIdentical(rep.int(NA, length(galp)), mcols(left(galp))[["FO"]])
}

test_readBamSortedPairs_obeyQname <- function()
{
    bf <- BamFile(fl, index=character(0), yieldSize=1, obeyQname=TRUE)
    scn1 <- scanBam(bf, obeyQname=TRUE)
    checkTrue(length(scn1[[1]]$qname) == 2)
    checkTrue(length(unique(scn1[[1]]$qname)) == 1)

    bf <- BamFile(fl, index=character(0), yieldSize=2, obeyQname=TRUE)
    scn2 <- scanBam(bf, obeyQname=TRUE)
    checkIdentical(scn1, scn2)

    bf <- BamFile(fl, index=character(0), yieldSize=3, obeyQname=TRUE)
    scn3 <- scanBam(bf, obeyQname=TRUE)
    checkTrue(length(scn3[[1]]$qname) == 4)
    checkTrue(length(unique(scn3[[1]]$qname)) == 2)
}
