library(pasillaBamSubset)
chr4 <- sortBam(untreated3_chr4(), tempfile(), byQname=TRUE)

test_readGAlignmentsListFromBam_noYieldSize <- function()
{
    fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
    fl_sort <- sortBam(fl, tempfile(), byQname=TRUE)
    bf <- BamFile(fl_sort, character(0), obeyQname=TRUE)
    galist <- readGAlignmentsListFromBam(fl_sort, group.as.pairs=FALSE)
    param <- ScanBamParam(what=c("rname", "strand", "pos", "cigar"),
                          flag=scanBamFlag(isUnmappedQuery=FALSE))
    scn0 <- scanBam(bf, param=param)
    checkTrue(validObject(galist))
    checkTrue(length(scn0[[1]]$rname) == sum(elementLengths(galist)))
}

test_readGAlignmentsListFromBam_yieldSize <- function()
{
    bf <- BamFile(chr4, character(0), yieldSize=1, obeyQname=TRUE)
    scn1 <- scanBam(bf)
    galist1 <- readGAlignmentsList(bf, group.as.pairs=FALSE)
    checkTrue(length(scn1[[1]]$qname) == 2)
    checkTrue(length(unique(scn1[[1]]$qname)) == 1)
    checkTrue(length(unique(scn1[[1]]$qname)) == length(galist1))

    bf <- BamFile(chr4, character(0), yieldSize=2, obeyQname=TRUE)
    scn2 <- scanBam(bf)
    galist2 <- readGAlignmentsList(bf, group.as.pairs=FALSE)
    checkTrue(length(scn2[[1]]$qname) == 4)
    checkTrue(length(unique(scn2[[1]]$qname)) == 2)
    checkTrue(length(unique(scn2[[1]]$qname)) == length(galist2))

    ## yieldSize > file length
    bf <- BamFile(chr4, character(0), yieldSize=176000, obeyQname=TRUE)
    scn3 <- scanBam(bf)
    checkTrue(length(scn3[[1]]$qname) == 175346)
    checkTrue(length(unique(scn3[[1]]$qname)) == 93620)

    ## yieldSize with 'while'
    bf <- open(BamFile(chr4, character(0), yieldSize=50000, obeyQname=TRUE))
    len <- 0 
    while (scn4 <- length(scanBam(bf)[[1]]$qname))
        len <- len + scn4
    close(bf)
    checkTrue(len == length(scn3[[1]]$qname))

    bf <- open(BamFile(chr4, character(0), yieldSize=50000, obeyQname=FALSE))
    len <- 0 
    while (scn5 <- length(scanBam(bf)[[1]]$qname))
        len <- len + scn5
    close(bf)
    checkTrue(len == countBam(bf)$records)
}

test_readGAlignmentsListFromBam_mcols <- function()
{
    bf <- BamFile(chr4, character(0), yieldSize=100, obeyQname=TRUE)
    param <- ScanBamParam(tag=("NM"))
    galist <- readGAlignmentsListFromBam(bf, param=param, group.as.pairs=FALSE)
    checkTrue(colnames(mcols(unlist(galist))) == "NM")
    param <- ScanBamParam(tag=("FO"))
    galist <- readGAlignmentsListFromBam(bf, param=param, group.as.pairs=FALSE)
    checkIdentical(rep.int(NA, length(unlist(galist))), 
                   mcols(unlist(galist))[["FO"]])
}

test_readGAlignmentsListFromBam_obeyQname <- function()
{
    bf <- BamFile(chr4, character(0), yieldSize=100, obeyQname=FALSE)
    checkException(readGAlignmentsListFromBam(bf), silent=TRUE)
}

test_readGAlignmentsListFromBam_groupAsPaired <- function()
{
    fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
    param <- ScanBamParam(what=c("flag", "mrnm", "mpos"))
    gal <- readGAlignmentsFromBam(fl, use.names=TRUE, param=param) 
    gal_paired <- groupAsPairs(gal)
    galist <- readGAlignmentsListFromBam(fl, param=param,
                                         group.as.pairs=FALSE)
    galist_paired <- groupAsPairs(galist)
    galp <- readGAlignmentPairsFromBam(fl, use.names=TRUE)

    ## groupAsPairs,GAlignments vs groupAsPairs,GAlignmentsList 
    checkTrue(length(gal_paired) == length(galist_paired))
    checkTrue(all(unique(names(gal_paired)) %in% unique(names(galist_paired))))

    ## groupAsPairs,GAlignments vs GAlignmentPairs
    ulist <- unlist(gal_paired, use.names=FALSE)
    paired <- ulist[mcols(ulist)$paired]
    checkTrue(all(names(paired) %in% names(galp)))
    checkTrue(length(paired) == 2*length(galp))
 
    ## groupAsPairs,GAlignmentsList vs GAlignmentPairs
    ulist <- unlist(galist_paired, use.names=FALSE)
    paired <- ulist[mcols(ulist)$paired]
    checkTrue(all(names(paired) %in% names(galp)))
    checkTrue(length(paired) == 2*length(galp))

    ## groupAsPairs,GAlignmentsList vs 'group.as.pairs=TRUE'
    galist2 <- readGAlignmentsListFromBam(fl, group.as.pairs=TRUE)
    ulist2 <- unlist(galist_paired, use.names=FALSE)
    paired2 <- ulist[mcols(ulist)$paired]
    ulist1 <- unlist(galist_paired, use.names=FALSE)
    paired1 <- ulist[mcols(ulist)$paired]
    checkTrue(all(names(paired1) %in% names(paired2)))
    checkTrue(length(paired1) == length(paired2))
}
