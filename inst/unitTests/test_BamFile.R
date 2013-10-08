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

test_BamFile_isIncomplete <- function()
{
    bf <- BamFile(fl, yieldSize=3000)
    checkIdentical(FALSE, isIncomplete(bf))
    open(bf)
    checkIdentical(TRUE, isIncomplete(bf))
    checkIdentical(3000L, length(scanBam(bf)[[1]][[1]]))
    checkIdentical(TRUE, isIncomplete(bf))
    checkIdentical(307L, length(scanBam(bf)[[1]][[1]]))
    checkIdentical(FALSE, isIncomplete(bf))
    close(bf)
    checkIdentical(FALSE, isIncomplete(bf))
}

test_BamFile_corrupt_index <- function()
{
    fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
    idx <- system.file("unitTests", "cases", "ex1_zero_index.bam.bai")
    bf <- BamFile(fl, idx)
    checkException(suppressWarnings(open(bf)), silent=TRUE)
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

test_BamFile_obeyQname <- function()
{
    srt <- sortBam(fl, tempfile(), byQname=TRUE)
    bf <- BamFile(srt, character(0), yieldSize=1, obeyQname=TRUE)
    scn <- scanBam(bf)
    checkTrue(length(scn[[1]]$qname) == 2)
    checkTrue(length(unique(scn[[1]]$qname)) == 1)

    bf <- BamFile(srt, character(0), yieldSize=2, obeyQname=TRUE)
    scn <- scanBam(bf)
    checkTrue(length(scn[[1]]$qname) == 4)
    checkTrue(length(unique(scn[[1]]$qname)) == 2)

    ## yieldSize is number of unique qnames
    ys <- 1000
    bf <- BamFile(srt, character(0), yieldSize=ys, obeyQname=TRUE)
    scn <- scanBam(bf)
    checkTrue(length(unique(scn[[1]]$qname)) == ys)

    ## yieldSize with 'while'
    bf <- open(BamFile(srt, character(0), yieldSize=500, obeyQname=TRUE))
    it <- integer() 
    while (length(res <- scanBam(bf)[[1]]$qname))
        it <- append(it, length(res))
    close(bf)
    checkIdentical(c(974L, 966L, 977L, 390L), it)
}

test_BamFile_asMates_all <- function()
{
    scn <- scanBam(fl)[[1]]
    scnm <- scanBam(BamFile(fl, asMates=TRUE))[[1]]
    param <- ScanBamParam()
    galp <- readGAlignmentPairsFromBam(fl, use.names=TRUE)

    ## qname
    checkTrue(all(scn$qname %in% scnm$qname))
    matenames <- scnm$qname[scnm$mates] 
    checkTrue(all(names(galp) %in% matenames))

    ## non-mates off last
    max1 <- max(which(scnm$mates))
    min0 <- min(which(!scnm$mates))
    checkTrue(max1 < min0)

    ## match GAlignmentPairs
    flag <- scanBamFlag(isPaired=TRUE, hasUnmappedMate=FALSE,
                        isUnmappedQuery=FALSE)
    param <- ScanBamParam(flag=flag, what=scanBamWhat())
    scnm <- scanBam(BamFile(fl, asMates=TRUE), param=param)[[1]]
    checkTrue(all(scnm$qname %in% names(galp)))
    checkTrue(all(names(galp) %in% scnm$qname))

    ## yieldSize - subset
    ## asMates=TRUE returns first 5 mates found
    ## asMates=FALSE returns first 5 records
    flag <- scanBamFlag(isPaired=TRUE, hasUnmappedMate=FALSE,
                        isUnmappedQuery=FALSE)
    param=ScanBamParam(flag=flag, what=scanBamWhat())
    scn <- scanBam(BamFile(fl, yieldSize=5), param=param)[[1]]
    scnm <- scanBam(BamFile(fl, asMates=TRUE, yieldSize=5), param=param)[[1]]
    checkTrue(length(scn$qname) == 5)
    checkTrue(length(scnm$qname) == 10)

    ## yieldSize - all records 
    scnm1 <- scanBam(BamFile(fl, asMates=TRUE))[[1]]
    scnm2 <- scanBam(BamFile(fl, asMates=TRUE, yieldSize=2000))[[1]]
    checkTrue(all(scnm1$qname %in% scnm2$qname))
    checkTrue(all(scnm1$mates == scnm2$mates))
    checkTrue(all(scnm1$.partition == scnm2$.partition))
}

test_BamFile_asMates_range <- function()
{
    which <- GRanges("seq1", IRanges(100, width=5))
    param <- ScanBamParam(which=which, what=scanBamWhat())
    scn <- scanBam(fl, param=param)[[1]]
    scnm <- scanBam(BamFile(fl, asMates=TRUE), param=param)[[1]]
    param <- ScanBamParam(which=which)
    galp <- readGAlignmentPairsFromBam(fl, param=param, use.names=TRUE)

    ## qname
    checkTrue(all(scn$qname %in% scnm$qname))
    matenames <- scnm$qname[scnm$mates] 
    checkTrue(all(names(galp) %in% matenames))
    checkTrue(length(scnm$qname) == sum(scnm$.partition))

    ## non-mates off last
    max1 <- max(which(scnm$mates))
    min0 <- min(which(!scnm$mates))
    checkTrue(max1 < min0)

    ## range - subset
    ## asMates=FALSE returns records in range
    ## asMates=TRUE mates all records in range 
    flag <- scanBamFlag(isPaired=TRUE, hasUnmappedMate=FALSE,
                        isUnmappedQuery=FALSE)
    param <- ScanBamParam(which=which, flag=flag, what=scanBamWhat())
    scn <- scanBam(fl, param=param)[[1]]
    scnm <- scanBam(BamFile(fl, asMates=TRUE), param=param)[[1]]
    checkTrue(length(scn$qname) == 8)
    checkTrue(length(scnm$qname) == 16)
    checkTrue(sum(scnm$mates) == 16)

    ## range - all records
    which <- GRanges(c("seq1", "seq2"), IRanges(start=c(0, 0), end=c(3000, 3000)))
    param <- ScanBamParam(which=which, what=scanBamWhat())
    scnm1 <- scanBam(BamFile(fl, asMates=TRUE))[[1]] ## get all w/ no param
    scnm2<- scanBam(BamFile(fl, asMates=TRUE), param=param) ## get all w/ param
    range2A <- scnm2[[1]]
    range2B <- scnm2[[2]]

    len_range <- length(range2A$qname) + length(range2B$qname)
    checkTrue(len_range == length(scnm1$qname))
    mates_range <- sum(range2A$mates, range2B$mates)
    checkTrue(mates_range == sum(scnm1$mates))
    mates_partition <- sum(range2A$.partition, range2B$.partition)
    checkTrue(mates_partition == sum(scnm1$.partition))
}

