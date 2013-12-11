fl <- system.file("extdata", "example.gtf.gz", package="Rsamtools")

test_TabixFile_openclose <- function()
{
    .normalizePath <- Rsamtools:::.normalizePath
    tab <- TabixFile(fl)
    checkTrue(!isOpen(tab))
    open(tab)
    checkTrue(isOpen(tab))

    checkIdentical(.normalizePath(fl), path(tab))
    checkIdentical(sprintf("%s.tbi", .normalizePath(fl)), index(tab))
    close(tab)
    checkTrue(!isOpen(tab))
    checkException(close(tab), silent=TRUE)
    tab <- open(tab)                   # open a closed TabixFile
    checkTrue(isOpen(tab))
    tab1 <- open(tab)                  # (re)open TabixFile
    checkTrue(isOpen(tab1))
    checkTrue(identical(tab$.extptr, tab1$.extptr))
}

test_TabixFile_scan <- function()
{
    tab <- open(TabixFile(fl))
    param <- GRanges("chr1", IRanges(1, 10000))

    res <- scanTabix(tab, param=param)
    checkIdentical("chr1:1-10000", names(res))
    checkIdentical(86L, length(res[[1]]))

    param <- GRanges(c("chr1", "chr2"),
                     IRanges(c(1,1), width=100000))

    res <- scanTabix(tab, param=param)
    checkIdentical(c("chr1:1-100000", "chr2:1-100000"), names(res))
    checkIdentical(c(157L, 15L), as.vector(sapply(res, length)))
}

test_TabixFile_yield <- function()
{
    tab <- open(TabixFile(fl, yieldSize=100))
    it <- integer()
    while(length(res <- scanTabix(tab)[[1]]))
        it <- append(it, length(res))
    close(tab)
    checkIdentical(c(100L, 100L, 37L), it)

    rng <- GRanges(c("seq1", "seq2"), IRanges(1, c(1575, 1584)))
    open(tab)
    checkException(scanTabix(tab, param=rng), silent=TRUE)
    close(tab)
}

test_TabixFile_header <- function()
{
    hdr <- headerTabix(fl)
    exp <- c("seqnames", "indexColumns", "skip", "comment", "header")
    checkIdentical(exp, names(hdr))
    checkIdentical(c("chr1", "chr2"), hdr$seqnames)
    checkIdentical(c(1L, 4L, 5L), unname(hdr$indexColumns))
    checkIdentical(0L, hdr$skip)
    checkIdentical("#", hdr$comment)
    checkIdentical(hdr$seqnames, seqnamesTabix(fl))
    checkIdentical(character(), hdr$header)
}

test_TabixFile_header_respects_first_line <- function()
{
    tbx <- open(TabixFile(fl))
    xx <- headerTabix(tbx)
    obs <- unlist(scanTabix(tbx))
    close(tbx)
    exp <- unlist(scanTabix(tbx))
    checkIdentical(obs, exp)
}

test_TabixFile_header_remote <- function()
{
    if ("windows" == .Platform$OS.type) {
        DEACTIVATED("remote tabix not supported on Windows")
        return(TRUE)
    }

    fl <- sprintf("%s/%s",
        "http://1000genomes.s3.amazonaws.com/release",
        "20110521/ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz")

    if (!is.null(tryCatch({open(con <- url(fl)); close(con)},
                 error=function(...) FALSE)))
        return(TRUE)

    obs <- headerTabix(fl)
    checkIdentical("22", obs$seqnames)
    exp <- structure(c(1L, 2L, 0L), .Names = c("seq", "start", "end"))
    checkIdentical(exp, obs$indexColumns)
    checkIdentical(0L, obs$skip)
    checkIdentical("#", obs$comment)
    checkIdentical(30L, length(obs$header))
}
