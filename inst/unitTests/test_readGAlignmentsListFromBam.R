library(pasillaBamSubset)
chr4 <- untreated3_chr4()

test_readGAlignmentsListFromBam_noYieldSize <- function()
{
    fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
    bf <- BamFile(fl, asMates=TRUE)
    galist <- readGAlignmentsListFromBam(fl)
    param <- ScanBamParam(what=c("rname", "strand", "pos", "cigar"),
                          flag=scanBamFlag(isUnmappedQuery=FALSE))
    scn0 <- scanBam(bf, param=param)
    checkTrue(validObject(galist))
    checkTrue(length(scn0[[1]]$rname) == sum(elementLengths(galist)))
}

test_readGAlignmentsListFromBam_yieldSize <- function()
{
    bf <- BamFile(chr4, asMates=TRUE, yieldSize=1)
    scn1 <- scanBam(bf)
    galist1 <- readGAlignmentsList(bf)
    checkTrue(length(scn1[[1]]$qname) == 2)
    checkTrue(length(unique(scn1[[1]]$qname)) == 1)
    checkTrue(length(unique(scn1[[1]]$qname)) == length(galist1))

    bf <- BamFile(chr4, asMates=TRUE, yieldSize=2)
    scn2 <- scanBam(bf)
    galist2 <- readGAlignmentsList(bf)
    checkTrue(length(scn2[[1]]$qname) == 4)
    checkTrue(length(unique(scn2[[1]]$qname)) == 2)
    checkTrue(length(unique(scn2[[1]]$qname)) == length(galist2))

    ## yieldSize > file length
    bf <- BamFile(chr4, asMates=TRUE, yieldSize=176000)
    scn3 <- scanBam(bf)
    checkTrue(length(scn3[[1]]$qname) == 175346)
    checkTrue(length(unique(scn3[[1]]$qname)) == 93620)

    ## yieldSize with 'while'
    bf <- open(BamFile(chr4, yieldSize=50000))
    len <- 0 
    while (scn4 <- length(scanBam(bf)[[1]]$qname))
        len <- len + scn4
    close(bf)
    checkTrue(len == length(scn3[[1]]$qname))

    bf <- open(BamFile(chr4, yieldSize=50000))
    len <- 0 
    while (scn5 <- length(scanBam(bf)[[1]]$qname))
        len <- len + scn5
    close(bf)
    checkTrue(len == countBam(bf)$records)
}

test_readGAlignmentsListFromBam_mcols <- function()
{
    bf <- BamFile(chr4, asMates=TRUE, yieldSize=100)
    param <- ScanBamParam(tag=("NM"))
    galist <- readGAlignmentsListFromBam(bf, param=param)
    current <- colnames(mcols(unlist(galist)))
    target <- c("mates", "NM")
    checkTrue(all(current %in% target))
    param <- ScanBamParam(tag=("FO"))
    galist <- readGAlignmentsListFromBam(bf, param=param)
    checkIdentical(rep.int(NA, length(unlist(galist))), 
                   mcols(unlist(galist))[["FO"]])
}

test_readGAlignmentsListFromBam_obeyQname <- function()
{
    bf <- BamFile(chr4, character(0), yieldSize=100, obeyQname=FALSE)
    galist <- readGAlignmentsListFromBam(bf)
    checkTrue(length(galist) == 100)
}

