library(Rsamtools)
library(RUnit)

fl <- system.file("extdata", "ex1.bam", package="Rsamtools")

.check0 <- function(res)
{
    checkIdentical("list", class(res))
    checkIdentical(scanBamWhat(), names(res))
    checkTrue(all(sapply(res, validObject)))
}

.check1 <- function(res)
{
    .check0(res)
    exp <- c("character", "integer", "factor", "factor", "integer",
             "integer", "integer", "Cigar", "factor", "integer",
             "integer", "DNAStringSet", "PhredQuality")
    checkIdentical(exp, as.vector(sapply(res, class)))
    checkIdentical(c("-", "+", "*"), levels(res[["strand"]]))
}

test_scanBam <- function()
{
    res <- scanBam(fl)[[1]]
    checkIdentical("list", class(res))
    checkIdentical(3307L, unique(sapply(res, length)))
    .check1(res)

    exp <- structure(c(11L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L), .Dim
                     = 10L, .Dimnames = structure(list( c("1", "2",
                     "3", "4", "6", "36", "37", "112", "283", "2804"
                     )), .Names = ""), class = "table")
                     checkIdentical(exp,
                     table(table(cigars(res[["cigar"]]))))

    exp <- structure(c(1641L, 1666L, 0L), .Dim = 3L, .Dimnames =
                     structure(list( c("-", "+", "*")), .Names = ""),
                     class = "table")
    checkIdentical(exp, table(res[["strand"]]))

    exp <- structure(c(8L, 40L, 858L, 17L, 714L, 5L, 12L, 11L, 35L,
                       714L, 18L, 858L, 12L, 5L), .Dim = 14L,
                       .Dimnames = structure(list(c("69", "73", "83",
                       "89", "99", "117", "121", "133", "137", "147",
                       "153", "163", "181", "185")), .Names = ""),
                       class = "table")
    checkIdentical(exp, table(res[["flag"]]))

    exp <- structure(c(6L, 37L, 2865L, 285L, 114L), .Dim = 5L,
                     .Dimnames = structure(list( c("33", "34", "35",
                     "36", "40")), .Names = ""), class = "table")
    checkIdentical(exp, table(width(res[["seq"]])))
    exp <- structure(c(36321L, 22060L, 22073L, 35958L, 139L), .Names =
                     c("A", "C", "G", "T", "other"))
    checkIdentical(exp,
                   alphabetFrequency(res[["seq"]], collapse=TRUE,
                                     baseOnly=TRUE))

    exp <- structure(c(263, 1, 20, 178, 577, 163, 195, 287, 286, 853,
                       290, 367, 424, 340, 395, 604, 601, 694, 898,
                       784, 1215, 2298, 1814, 1889, 3330, 10633,
                       80537, 6444, 150, 18, 3), .Names = c("!", "#",
                       "$", "%", "&", "'", "(", ")", "*", "+", ",",
                       "-", ".", "/", "0", "1", "2", "3", "4", "5",
                       "6", "7", "8", "9", ":", ";", "<", "=", ">",
                       "?", "@"))
    checkIdentical(exp, rowSums(consensusMatrix(res[["qual"]])))
}

test_scanBam_which<- function()
{
    ## 'which'
    which <- RangesList(seq1=IRanges(1000, 2000),
                    seq2=IRanges(c(100, 1000), c(1000, 2000)))
    p1 <- ScanBamParam(which=which)
    res <- scanBam(fl, param=p1)

    checkIdentical("list", class(res))
    exp <- c("seq1:1000-2000", "seq2:100-1000", "seq2:1000-2000")
    checkIdentical(exp, names(res))
    for (i in seq_along(res)) .check1(res[[i]])

    exp <- structure(c(612L, 1168L, 642L),
                     .Names = c("seq1:1000-2000", "seq2:100-1000",
                     "seq2:1000-2000"))
    checkIdentical(exp,
                   sapply(res, function(x) unique(sapply(x, length))))

}

test_scanBam_what <- function()
{
    p2 <- ScanBamParam(what=c("rname", "strand", "pos", "width"))
    res <- scanBam(fl, param=p2)
    checkIdentical("list", class(res))
    checkIdentical(1L, length(res))
    res1 <- res[[1]]
    .check0(res1)

    checkIdentical(9L, sum(sapply(res1, is.null)))

    exp <- structure(c("factor", "factor", "integer", "integer"),
                     .Names = c("rname", "strand", "pos", "width"))
    checkIdentical(exp, sapply(res1[3:6], class))
    checkIdentical(3307L, unique(sapply(res1[3:6], length)))
}

test_scanBam_flag <- function()
{
    p3 <- ScanBamParam(flag=scanBamFlag(isMinusStrand=TRUE))
    res <- scanBam(fl, param=p3)
    checkIdentical("list", class(res))
    checkIdentical(1L, length(res))
    res1 <- res[[1]]
    .check1(res1)

    checkIdentical(1641L, unique(sapply(res1, length)))
}

test_scanBam()
test_scanBam_which()
test_scanBam_what()
test_scanBam_flag()


