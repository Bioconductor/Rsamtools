fl <- system.file("extdata", "ex1.bcf", package="Rsamtools")

test_BcfFile_openclose <- function()
{
    .normalizePath <- Rsamtools:::.normalizePath
    bf <- BcfFile(fl, character(0))     # no index
    checkTrue(!isOpen(bf))
    open(bf)
    checkTrue(isOpen(bf))
    checkIdentical(.normalizePath(fl), path(bf))
    checkIdentical(character(0), index(bf))
    close(bf)
    checkTrue(!isOpen(bf))
    checkException(close(bf), silent=TRUE)
    bf <- open(bf)                   # open a closed BcfFile
    checkTrue(isOpen(bf))
    bf1 <- open(bf)                  # (re)open BcfFile
    checkTrue(isOpen(bf1))
    checkTrue(identical(bf$.extptr, bf1$.extptr))

    ## checkTrue(isOpen(bf))
    ## checkIdentical(.normalizePath(fl), path(bf))
    ## checkIdentical(.normalizePath(fl), index(bf))
}

test_BcfFile_scanBcfHeader <- function()
{
    .chk <- function(h) {
        checkTrue(validObject(h))
        checkEquals(3L, length(h))
        checkEquals(2L, length(h[["Reference"]]))
        checkEquals("ex1.bam", h[["Sample"]])
        checkEquals(1L, length(h[["Header"]]))
        checkEquals(c(0, 0), dim(h[["Header"]][["META"]]))
    }
    bf <- open(BcfFile(fl, character(0)))
    .chk(scanBcfHeader(bf))
    bf <- BcfFile(fl, character(0))
    .chk(scanBcfHeader(bf))
    checkTrue(!isOpen(bf))

    header <-
        structure(list(Reference = character(0), Sample =
                       character(0), Header =
                       "##FORMAT=<ID=GD,Number=1,Type=Float,Description=\"alleles [0,2]\">"),
                       .Names = c("Reference", "Sample", "Header"))
    checkTrue(validObject(Rsamtools:::.bcfHeaderAsSimpleList(header)))
}

test_BcfFile_scan_noindex <- function()
{
    bf <- open(BcfFile(fl, character(0)))
    checkTrue(isOpen(bf))

    p <- ScanBcfParam()
    res <- scanBcf(bf, param=p)
    checkEquals(11L, length(res))
    checkEquals(3065L, res[["RecordsPerRange"]])
    checkEquals(3065L, unique(sapply(res[1:9], length)))
    checkEquals(1, length(res[["GENO"]]))
    checkEquals(3065L, length(res[["GENO"]][["PL"]]))

    checkIdentical(res, scanBcf(bf))
}

test_BcfFile_scan_index <- function()
{
    bf <- open(BcfFile(fl))
    which <- RangesList(seq1=IRanges(c(1, 1000), width=10),
                        seq2=IRanges(c(100, 1000), width=10))
    p <- ScanBcfParam(which=which)
    res <- scanBcf(bf, param=p)
    checkEquals(11L, length(res))
    checkEquals(c(0, 10, 10, 20), res[["RecordsPerRange"]])
    checkEquals(30L, unique(sapply(res[1:9], length)))
    checkEquals(1, length(res[["GENO"]]))
    checkEquals(30L, length(res[["GENO"]][["PL"]]))
}
