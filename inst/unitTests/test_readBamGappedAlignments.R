test_readBamGappedAlignments <- function()
{
    fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
    which <- RangesList(seq1=IRanges(1, width=100))
    param <- ScanBamParam(which=which)
    result <- readBamGappedAlignments(fl, param=param)
    checkTrue(validObject(result))
    checkIdentical(c(seq1=1575L, seq2=1584L), seqlengths(result))
}

test_readBamGappedAlignments_missing_param <- function()
{
    fl <- system.file("unitTests", "cases", "ex1_noindex.bam",
                      package="Rsamtools")
    result0 <- readBamGappedAlignments(fl)
    checkTrue(validObject(result0))

    bf <- open(BamFile(fl, character()))
    result1 <- readBamGappedAlignments(bf)
    checkIdentical(result1, result0)
}

test_readBamGappedAlignments_length0 <- function()
{
    fl <- system.file("extdata", "ex1.bam", package="Rsamtools")

    which <- RangesList(seq1=IRanges(100000, width=100))
    param <- ScanBamParam(which=which)
    result <- readBamGappedAlignments(fl, param=param)
    checkTrue(validObject(result))

    which <- RangesList(seq1=IRanges(c(1, 100000), width=100))
    param <- ScanBamParam(which=which)
    result <- readBamGappedAlignments(fl, param=param)
    checkTrue(validObject(result))
}

test_readBamGappedAlignments_tag <- function()
{
    fl <- system.file("extdata", "ex1.bam", package="Rsamtools")

    ## valid
    param <- ScanBamParam(tag=("NM"))
    aln <- readBamGappedAlignments(fl, param=param)
    checkIdentical(924L, sum(mcols(aln)[["NM"]]))

    ## empty
    param <- ScanBamParam(tag=("FO"))
    aln <- readBamGappedAlignments(fl, param=param)
    checkIdentical(rep.int(NA, length(aln)), mcols(aln)[["FO"]])
}
