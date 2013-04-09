test_readGAlignmentsFromBam <- function()
{
    fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
    which <- RangesList(seq1=IRanges(1, width=100))
    param <- ScanBamParam(which=which)
    result <- readGAlignmentsFromBam(fl, param=param)
    checkTrue(validObject(result))
    checkIdentical(c(seq1=1575L, seq2=1584L), seqlengths(result))
}

test_readGAlignmentsFromBam_missing_param <- function()
{
    fl <- system.file("unitTests", "cases", "ex1_noindex.bam",
                      package="Rsamtools")
    result0 <- readGAlignmentsFromBam(fl)
    checkTrue(validObject(result0))

    bf <- open(BamFile(fl, character()))
    result1 <- readGAlignmentsFromBam(bf)
    checkIdentical(result1, result0)
}

test_readGAlignmentsFromBam_length0 <- function()
{
    fl <- system.file("extdata", "ex1.bam", package="Rsamtools")

    which <- RangesList(seq1=IRanges(100000, width=100))
    param <- ScanBamParam(which=which)
    result <- readGAlignmentsFromBam(fl, param=param)
    checkTrue(validObject(result))

    which <- RangesList(seq1=IRanges(c(1, 100000), width=100))
    param <- ScanBamParam(which=which)
    result <- readGAlignmentsFromBam(fl, param=param)
    checkTrue(validObject(result))
}

test_readGAlignmentsFromBam_tag <- function()
{
    fl <- system.file("extdata", "ex1.bam", package="Rsamtools")

    ## valid
    param <- ScanBamParam(tag=("NM"))
    gal <- readGAlignmentsFromBam(fl, param=param)
    checkIdentical(924L, sum(mcols(gal)[["NM"]]))

    ## empty
    param <- ScanBamParam(tag=("FO"))
    gal <- readGAlignmentsFromBam(fl, param=param)
    checkIdentical(rep.int(NA, length(gal)), mcols(gal)[["FO"]])
}
