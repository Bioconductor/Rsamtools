library(pasillaBamSubset)
fl <- sortBam(untreated3_chr4(), tempfile(), byQname=TRUE)

test_readBamGAlignments_pairs <- function()
{
    bf <- BamFile(fl, index=character(0), yieldSize=100, obeyQname=TRUE)
    galp <- readBamGappedAlignmentPairs(bf)
    checkTrue(validObject(galp))
    checkTrue(length(galp) == 50)
}

test_readBamGAlignmentsList_length <- function()
{
    bf <- BamFile(fl, index=character(0), yieldSize=100, obeyQname=TRUE)
    galist <- readBamGAlignmentsList(bf, asProperPairs=FALSE)
    checkTrue(validObject(galist))
    checkTrue(length(galist) == 51)

    bf <- BamFile(fl, index=character(0), yieldSize=3, obeyQname=TRUE)
    galist <- readBamGAlignmentsList(bf, asProperPairs=FALSE)
    checkTrue(validObject(galist))
    checkTrue(length(galist) == 2)
}

test_readBamGAlignmentsList_mcols <- function()
{
    bf <- BamFile(fl, index=character(0), yieldSize=100, obeyQname=TRUE)
    param <- ScanBamParam(tag=("NM"))
    galist <- readBamGAlignmentsList(bf, param=param, asProperPairs=FALSE)
    checkTrue(colnames(mcols(unlist(galist))) == "NM")
    param <- ScanBamParam(tag=("FO"))
    galist <- readBamGAlignmentsList(bf, param=param, asProperPairs=FALSE)
    checkIdentical(rep.int(NA, length(unlist(galist))), 
                   mcols(unlist(galist))[["FO"]])
}

test_readBamGAlignmentsList_obeyQname <- function()
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
