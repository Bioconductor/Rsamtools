## suppressMessages({
##     library(Rsamtools)
##     library(RUnit)
## })

fl <- system.file(package="Rsamtools", "extdata", "tagfilter.bam")
bf <- BamFile(fl)

msbp <- function(rnames) {
    ScanBamParam(which=GRanges(rnames, IRanges(1,10)), what="rname")
}
numrecs <- function(sbp) countBam(fl, param=sbp)[["records"]]

test_notags <- function() {
    sbp <- msbp("notags")
    bamTagFilter(sbp) <- list(TT="bogus")
    res <- numrecs(sbp)
    checkIdentical(0L, res)
}
##test_notags()

test_shared_by_multiple_reads <- function() {
    sbp <- ScanBamParam(tagFilter=list(AA=c("a", "d")))
    res <- numrecs(sbp)
    checkIdentical(2L, res)
}
##test_shared_by_multiple_reads()

test_multitags <- function() {
    ## two reads have AA:A:a or AA:A:d, but only one of the two has
    ## II:i:45
    sbp <- ScanBamParam(tagFilter=list(AA=c("a", "d"), II=45))
    res <- numrecs(sbp)
    checkIdentical(1L, res)
}
##test_multitags()

test_integer <- function() {
    sbp <- msbp("itag")
    ## exclude all
    bamTagFilter(sbp) <- list(II=1)
    res <- numrecs(sbp)
    checkIdentical(0L, res)

    ## include 2 discontiguous
    bamTagFilter(sbp) <- list(II=c(42, 44))
    res <- numrecs(sbp)
    checkIdentical(2L, res)

    ## exception for mismatch
    bamTagFilter(sbp) <- list(II="fun")
    checkException(numrecs(sbp))
}
##test_integer()

## per the SAM spec, single printable character is different from a
## string
test_single_printable <- function() {
    sbp <- msbp("Atag")
    ## exclude all
    bamTagFilter(sbp) <- list(AA="d")
    res <- numrecs(sbp)
    checkIdentical(0L, res)

    ## include 2 discontiguous
    bamTagFilter(sbp) <- list(AA=c("a", "c"))
    res <- numrecs(sbp)
    checkIdentical(2L, res)

    ## exception for mismatch
    bamTagFilter(sbp) <- list(AA="fun")
    checkException(numrecs(sbp))
}
##test_single_printable()

test_string <- function() {
    sbp <- msbp("Ztag")
    ## exclude all
    bamTagFilter(sbp) <- list(ZZ="wok")
    res <- numrecs(sbp)
    checkIdentical(0L, res)

    ## include 2 discontiguous
    bamTagFilter(sbp) <- list(ZZ=c("woo", "wow"))
    res <- numrecs(sbp)
    checkIdentical(2L, res)

    ## exception for mismatch
    bamTagFilter(sbp) <- list(ZZ=1)
    checkException(numrecs(sbp))
}
##test_string()

## confirm throwing error when user tries to filter on a tag that has
## an unsupported type in the BAM file
test_unsupported_tag_types <- function() {
    ## floating point type
    checkException(countBam(fl, param=ScanBamParam(tagFilter=list(FF=13))))
    ## hex array
    checkException(countBam(fl, param=ScanBamParam(tagFilter=list(HH="foo"))))
    ## integer or numeric *array*
    checkException(countBam(fl, param=ScanBamParam(tagFilter=list(BB="foo"))))
}
##test_unsupported_tag_types()

## Input validation

test_exception_names <- function() {
    ## Too many letters
    checkException(ScanBamParam(tagFilter=list(NNN=1)))
    ## Too few
    checkException(ScanBamParam(tagFilter=list(N=1)))
    ## No names
    checkException(ScanBamParam(tagFilter=list(1)))
}
##test_exception_names()

test_exception_floating_point <- function() {
    checkException(ScanBamParam(tagFilter=list(FF=13.001)))
}
##test_exception_floating_point()

test_exception_weird_values <- function() {
    checkException(ScanBamParam(tagFilter=list(FF=NULL)))
    checkException(ScanBamParam(tagFilter=list(FF=NA)))
    ## zero-length
    checkException(ScanBamParam(tagFilter=list(FF=character())))
    ## empty string
    checkException(ScanBamParam(tagFilter=list(FF="")))
}
##test_exception_weird_values()
