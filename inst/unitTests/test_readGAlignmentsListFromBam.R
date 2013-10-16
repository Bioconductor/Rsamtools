library(pasillaBamSubset)
chr4 <- untreated3_chr4()

test_readGAlignmentsListFromBam_construction <- function()
{
    fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
    bf <- BamFile(fl, asMates=TRUE)
    galist <- readGAlignmentsListFromBam(fl)
    checkTrue(is.null(names(galist)))
    galist <- readGAlignmentsListFromBam(fl, use.names=TRUE)
    target <- c("EAS54_61:4:143:69:578", "EAS219_FC30151:7:51:1429:1043")
    checkIdentical(names(galist)[1:2], target)

    ## first segment first
    param <- ScanBamParam(what="flag")
    galist <- readGAlignmentsListFromBam(fl, param=param)
    mates <- galist[unlist(mcols(galist)$mates)]
    flagBit <- bamFlagAsBitMatrix(mcols(unlist(mates))$flag,
                                  bitnames="isFirstMateRead") 
    m <- matrix(flagBit, nrow=2)
    checkTrue(rowSums(m)[2] == 0L)
}

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
    checkTrue(names(mcols(galist)) == "mates")

    param <- ScanBamParam(tag=("FO"))
    galist <- readGAlignmentsListFromBam(bf, param=param)
    checkIdentical(rep.int(NA, length(unlist(galist))), 
                   mcols(unlist(galist))[["FO"]])
}
