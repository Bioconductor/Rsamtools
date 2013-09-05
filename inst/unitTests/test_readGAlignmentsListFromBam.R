library(pasillaBamSubset)
chr4 <- untreated3_chr4()

test_readGAlignmentsListFromBam_noYieldSize <- function()
{
    fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
    bf <- BamFile(fl, asMates=TRUE)
    galist <- readGAlignmentsListFromBam(fl)
    checkTrue(validObject(galist))
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
