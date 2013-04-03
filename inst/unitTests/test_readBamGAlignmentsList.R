library(pasillaBamSubset)
fl <- sortBam(untreated3_chr4(), tempfile(), byQname=TRUE)

test_readBamGAlignments_noYieldSize <- function()
{
    bf <- BamFile(fl, character(0), obeyQname=TRUE)
    galist <- readBamGAlignmentsList(bf, asProperPairs=FALSE)
    scn0 <- scanBam(bf)
    checkTrue(validObject(galist))
    checkTrue(length(scn0[[1]]$qname) == 
              length(unlist(galist, use.names=FALSE)))
}

test_readBamGAlignmentsList_yieldSize <- function()
{
    bf <- BamFile(fl, character(0), yieldSize=1, obeyQname=TRUE)
    scn1 <- scanBam(bf)
    galist1 <- readGAlignmentsList(bf, asProperPairs=FALSE)
    checkTrue(length(scn1[[1]]$qname) == 2)
    checkTrue(length(unique(scn1[[1]]$qname)) == 1)
    checkTrue(length(unique(scn1[[1]]$qname)) == length(galist1))

    bf <- BamFile(fl, character(0), yieldSize=2, obeyQname=TRUE)
    scn2 <- scanBam(bf)
    galist2 <- readGAlignmentsList(bf, asProperPairs=FALSE)
    checkTrue(length(scn2[[1]]$qname) == 4)
    checkTrue(length(unique(scn2[[1]]$qname)) == 2)
    checkTrue(length(unique(scn2[[1]]$qname)) == length(galist2))

    ## yieldSize > file length
    bf <- BamFile(fl, character(0), yieldSize=176000, obeyQname=TRUE)
    scn3 <- scanBam(bf)
    checkTrue(length(scn3[[1]]$qname) == 175346)
    checkTrue(length(unique(scn3[[1]]$qname)) == 93620)

    ## yieldSize with 'while'
    bf <- open(BamFile(fl, character(0), yieldSize=50000, obeyQname=TRUE))
    len <- 0 
    while (scn4 <- length(scanBam(bf)[[1]]$qname))
        len <- len + scn4
    close(bf)
    checkTrue(len == length(scn3[[1]]$qname))

    bf <- open(BamFile(fl, character(0), yieldSize=50000, obeyQname=FALSE))
    len <- 0 
    while (scn5 <- length(scanBam(bf)[[1]]$qname))
        len <- len + scn5
    close(bf)
    checkTrue(len == countBam(bf)$records)
}

test_readBamGAlignmentsList_mcols <- function()
{
    bf <- BamFile(fl, character(0), yieldSize=100, obeyQname=TRUE)
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
    ## TBD
}
