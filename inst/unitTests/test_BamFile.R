fl <- system.file("extdata", "ex1.bam", package="Rsamtools")

test_BamFile_openclose <- function()
{
    .normalizePath <- Rsamtools:::.normalizePath
    bf <- open(BamFile(fl))
    checkTrue(isOpen(bf))
    checkIdentical(.normalizePath(fl), path(bf))
    checkIdentical(sprintf("%s.bai", .normalizePath(fl)), index(bf))
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
    exp <- c(1000L, 1000L, 1000L, 307L)
    checkIdentical(exp, it)

    bf <- open(BamFile(fl, yieldSize=500))
    which <- GenomicRanges::tileGenome(seqlengths(bf), tilewidth=500,
                                       cut.last.tile.in.chrom=TRUE)
    param <- ScanBamParam(what="rname", which=which)
    it <- integer()
    while(res <- sum(unlist(scanBam(bf, param=param), use.names=FALSE)))
        it <- append(it, res)
    close(bf)
    checkIdentical(6L, length(it))
    ## over-count, since reads overlap more than one range
    checkIdentical(5472L, sum(it))
}

test_BamFileList_constructor <- function()
{
    checkTrue(validObject(res <- BamFileList(fl)))
    checkIdentical(fl, path(res[[1]]))
    checkIdentical(sprintf("%s.bai", fl), index(res[[1]]))
    
    checkTrue(validObject(res <- BamFileList(fl, fl)))
    checkIdentical(fl, path(res[[1]]))
    checkIdentical(fl, index(res[[1]]))
    ## old-style index=character()
    checkTrue(validObject(res <- BamFileList(fl, index=character())))
    checkIdentical(fl, path(res[[1]]))
    checkIdentical(character(), index(res[[1]]))

    ## vector of files
    fl <- c(fl, fl)
    checkTrue(validObject(res <- BamFileList(fl)))
    checkIdentical(fl, unname(path(res)))
    checkIdentical(sprintf("%s.bai", fl), unname(sapply(res, index)))

    checkTrue(validObject(res <- BamFileList(fl, fl)))
    checkIdentical(fl, unname(path(res)))
    checkIdentical(fl, unname(sapply(res, index)))
    ## old-style index=character()
    checkTrue(validObject(res <- BamFileList(fl, index=character())))
    checkIdentical(fl, unname(path(res)))
    checkIdentical(list(character(), character()), unname(lapply(res, index)))
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
    library(GenomicAlignments)
    param <- ScanBamParam()
    galp <- readGAlignmentPairsFromBam(fl, use.names=TRUE)

    ## qname
    checkTrue(all(scn$qname %in% scnm$qname))
    matenames <- scnm$qname[scnm$mate_status == "mated"] 
    checkTrue(all(names(galp) %in% matenames))

    ## non-mates off last
    max1 <- max(which(scnm$mate_status == "mated"))
    min0 <- min(which(scnm$mate_status != "mated"))
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
    checkTrue(all(scnm1$mate_status == scnm2$mate_status))
    checkTrue(all(scnm1$groupid == scnm2$groupid))
}

test_BamFile_asMates_range <- function()
{
    which <- GRanges("seq1", IRanges(100, width=5))
    param <- ScanBamParam(which=which, what=scanBamWhat())
    scn <- scanBam(fl, param=param)[[1]]
    scnm <- scanBam(BamFile(fl, asMates=TRUE), param=param)[[1]]
    library(GenomicAlignments)
    param <- ScanBamParam(which=which)
    galp <- readGAlignmentPairsFromBam(fl, param=param, use.names=TRUE)

    ## qname
    checkTrue(all(scn$qname %in% scnm$qname))
    matenames <- scnm$qname[scnm$mate_status == "mated"] 
    checkTrue(all(names(galp) %in% matenames))
    checkTrue(length(scnm$qname) == length(scnm$groupid))

    ## non-mates off last
    max1 <- max(which(scnm$mate_status == "mated"))
    min0 <- min(which(scnm$mate_status != "mated"))
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
    checkTrue(sum(scnm$mate_status == "mated") == 16)

    ## range - all records
    which <- GRanges(c("seq1", "seq2"), IRanges(start=c(0, 0), end=c(3000, 3000)))
    param <- ScanBamParam(which=which, what=scanBamWhat())
    scnm1 <- scanBam(BamFile(fl, asMates=TRUE))[[1]] ## get all w/ no param
    scnm2<- scanBam(BamFile(fl, asMates=TRUE), param=param) ## get all w/ param
    range2A <- scnm2[[1]]
    range2B <- scnm2[[2]]

    len_range <- length(range2A$qname) + length(range2B$qname)
    checkTrue(len_range == length(scnm1$qname))
    mates_range <- sum(range2A$mate_status == "mated", 
                       range2B$mate_status == "mated")
    checkTrue(mates_range == sum(scnm1$mate_status == "mated"))
    mates_partition <- length(range2A$.partition) + length(range2B$.partition)
    checkTrue(mates_partition == length(scnm1$.partition))
}

