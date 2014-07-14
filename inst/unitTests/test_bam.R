fl <- system.file("extdata", "ex1.bam", package="Rsamtools")

.check0 <- function(res)
{
    checkIdentical("list", class(res))
    checkTrue(all(names(res) %in% scanBamWhat()))
    checkTrue(all(sapply(res, validObject)))
}

.check1 <- function(res)
{
    .check0(res)
    exp <- c("character", "integer", "factor", "factor", "integer",
             "integer", "integer", "character", "factor", "integer",
             "integer", "DNAStringSet", "PhredQuality")
    checkIdentical(exp, as.vector(sapply(res, class)))
    checkIdentical(levels(strand()), levels(res[["strand"]]))
}

test_strand <- function()
{

    exp <- structure(integer(0), .Label = c("+", "-", "*"), class =
                     "factor")
    checkIdentical(exp, strand())
    checkIdentical(levels(exp), Rsamtools:::.STRAND_LEVELS)
}

test_scanBam <- function()
{
    res <- scanBam(fl)[[1]]
    checkIdentical("list", class(res))
    checkIdentical(3307L, unique(sapply(res, length)))
    .check1(res)

    exp <- structure(c(11L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L), .Dim
                     = 9L, .Dimnames = structure(list( c("1", "2",
                     "3", "4", "6", "37", "112", "283", "2804"
                     )), .Names = ""), class = "table")
    checkIdentical(exp, table(table(res[["cigar"]])))

    exp <- structure(c(1647L, 1624L, 0L, 36L), .Dim = 4L, .Dimnames =
                     structure(list( c(levels(strand()), NA)), .Names =
                     ""), class = "table")
    checkIdentical(exp, table(res[["strand"]], useNA="always"))

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
    exp <- structure(c(40144L, 23398L, 20735L, 32135L, 139L),
                     .Names = c("A", "C", "G", "T", "other"))
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
    checkIdentical(exp,
                   rowSums(consensusMatrix(res[["qual"]])))
}

test_scanBam_which <- function()
{
    ## 'which'
    which <- RangesList(seq1=IRanges(1000, 2000),
                    seq2=IRanges(c(100, 1000), c(1000, 2000)))
    p1 <- ScanBamParam(which=which, what=scanBamWhat())
    res <- scanBam(fl, param=p1)

    checkIdentical("list", class(res))
    exp <- c("seq1:1000-2000", "seq2:100-1000", "seq2:1000-2000")
    checkIdentical(exp, names(res))
    for (i in seq_along(res)) .check1(res[[i]])

    exp <- structure(c(612L, 1169L, 642L),
                     .Names = c("seq1:1000-2000", "seq2:100-1000",
                     "seq2:1000-2000"))
    checkIdentical(exp,
                   sapply(res, function(x) unique(sapply(x, length))))

}

test_scanBam_which_bounds <- function()
{
    ## 'which'
    which <- RangesList(seq1=IRanges(1000, 2000),
                    seq2=IRanges(c(100, 1000), c(1000, 2000)))
    p1 <- ScanBamParam(which=which, what=scanBamWhat())
    res <- scanBam(fl, param=p1)

    checkIdentical("list", class(res))
    exp <- c("seq1:1000-2000", "seq2:100-1000", "seq2:1000-2000")
    checkIdentical(exp, names(res))
    for (i in seq_along(res)) .check1(res[[i]])

    exp <- structure(c(612L, 1169L, 642L),
                     .Names = c("seq1:1000-2000", "seq2:100-1000",
                     "seq2:1000-2000"))
    checkIdentical(exp,
                   sapply(res, function(x) unique(sapply(x, length))))

    ## ranges correct?
    for (i in c(1, 2, 10, 30, 34, 35, 36, 37, 460)) {
        snp <- IRanges(i, i)
        which <- GRanges("seq2", IRanges(1, 1000))
        param <- ScanBamParam(which=which, what=c("pos", "qwidth"))
        res <- with(scanBam(fl, param=param)[[1]], {
            idx <- !is.na(pos)
            IRanges(pos[idx], width=qwidth[idx])
        })
        exp <- sum(countOverlaps(res, snp))
        param <- ScanBamParam(which=GRanges("seq2", snp), what=scanBamWhat())
        try(checkIdentical(exp, countBam(fl, param=param)$records,
                           msg=sprintf("i: %d\n", i)))
    }
}

test_scanBam_which_order <- function()
{
    which <- RangesList(seq2=IRanges(c(1000, 100), c(2000, 1000)),
                        seq1=IRanges(1000, 2000))
    p1 <- ScanBamParam(which=which, what="pos")
    res <- scanBam(fl, param=p1)
    checkIdentical("list", class(res))
    exp <- structure(c(642L, 1169L, 612L),
                     .Names = c( "seq2:1000-2000", "seq2:100-1000",
                       "seq1:1000-2000"))
    obs <- sapply(res, function(x) unique(sapply(x, length)))
    checkIdentical(exp, obs)
}

test_scanBam_which_empty <- function()
{
    ## range 1 is empty
    which <- IRangesList(seq2=IRanges(c(1570,1562), width=1))
    what <- c("strand", "rname", "mrnm")
    res <- scanBam(fl, param=ScanBamParam(what=what, which=which))

    checkIdentical(c(3L, 3L), unname(sapply(res, length)))
    checkIdentical(res[[1]][["strand"]], strand())
    exp <- factor(levels=c("seq1", "seq2"))
    checkIdentical(exp, res[[1]][["rname"]])
    checkIdentical(exp, res[[1]][["mrnm"]])

    checkIdentical(res[[2]][["strand"]], strand(rep("-", 3)))
    exp <- factor(rep("seq2", 3), levels=c("seq1", "seq2"))
    checkIdentical(exp, res[[2]][["rname"]])
    checkIdentical(exp, res[[2]][["mrnm"]])
}

test_scanBam_what_overflow <- function()
{
    ## src/samtools/bam_index.c requires that the largest bin be <= 512Mbp
    which <- RangesList(seq1=IRanges(1000, 536870912L-1L))
    p1 <- ScanBamParam(which=which, what=scanBamWhat())
    checkTrue(validObject(scanBam(fl, param=p1)))
    bamWhich(p1) <- RangesList(seq1=IRanges(1000, 536870912L))
    checkTrue(validObject(scanBam(fl, param=p1)))
    bamWhich(p1) <- RangesList(seq1=IRanges(1000, 536870912L+1L))
    xx <- tryCatch(scanBam(fl, param=p1), error=conditionMessage)
    checkIdentical("'end' must be <= 536870912", strsplit(xx, "\n")[[1]][1])
    checkTrue(validObject(scanBam(fl)));
}

test_scanBam_what <- function()
{
    p2 <- ScanBamParam(what=c("rname", "strand", "pos", "qwidth"))
    res <- scanBam(fl, param=p2)
    checkIdentical("list", class(res))
    checkIdentical(1L, length(res))
    res1 <- res[[1]]
    .check0(res1)

    checkIdentical(4L, length(res1))

    exp <- structure(c("factor", "factor", "integer", "integer"),
                     .Names = c("rname", "strand", "pos", "qwidth"))
    checkIdentical(exp, sapply(res1, class))
    checkIdentical(3307L, unique(sapply(res1, length)))
}

test_scanBam_flag <- function()
{
    p3 <- ScanBamParam(flag=scanBamFlag(isMinusStrand=TRUE),
                       what=scanBamWhat())
    res <- scanBam(fl, param=p3)
    checkIdentical("list", class(res))
    checkIdentical(1L, length(res))
    res1 <- res[[1]]
    .check1(res1)

    checkIdentical(1641L, unique(sapply(res1, length)))
}

test_scanBam_badSpace <- function()
{
    which <- RangesList(badspc=IRanges(100000, 2000000))
    p1 <- ScanBamParam(which=which, what=scanBamWhat())

    oopts <- options(warn=-1)
    on.exit(options(oopts))
    flag <- list()
    tryCatch({
        withCallingHandlers(scanBam(fl, param=p1), warning=function(w) {
            checkTrue(grepl("space 'badspc'", conditionMessage(w)))
            flag[["warn"]] <<- TRUE
        })
    }, error=function(e) {
        fl0 <- Rsamtools:::.normalizePath(fl)
        checkTrue(grepl(paste("file:", fl0), conditionMessage(e), fixed=TRUE))
        flag[["err"]] <<- TRUE
    }, finally=local({
        checkEquals(2L, length(flag))
        checkTrue(all(flag))
        flag[["flag"]] <<- TRUE
    }))
    checkTrue(flag[["flag"]])
}

test_scanBam_index <- function()
{
    which <- RangesList(seq1=IRanges(1000, 2000),
                    seq2=IRanges(c(100, 1000), c(1000, 2000)))
    p1 <- ScanBamParam(which=which, what=scanBamWhat())
    src <- system.file("unitTests", "cases", package="Rsamtools")
    fl <- file.path(src, "ex1_noindex.bam")
    idx <- system.file("extdata", "ex1.bam", package="Rsamtools")

    checkException(scanBam(fl, character(), param=p1), silent=TRUE)

    res <- scanBam(fl, idx, param=p1)
    checkIdentical("list", class(res))
    exp <- c("seq1:1000-2000", "seq2:100-1000", "seq2:1000-2000")
    checkIdentical(exp, names(res))
    for (i in seq_along(res)) .check1(res[[i]])

    exp <- structure(c(612L, 1169L, 642L),
                     .Names = c("seq1:1000-2000", "seq2:100-1000",
                     "seq2:1000-2000"))
    checkIdentical(exp,
                   sapply(res, function(x) unique(sapply(x, length))))

    checkException(suppressWarnings(scanBam(fl, tempfile(), 
                   param=p1), silent=TRUE))
}

test_scanBam_sam <- function()
{
    tbl <- read.table(system.file("extdata", "ex1.sam",
                                  package="Rsamtools"),
                      sep="\t", quote="", fill=TRUE, comment="")
    bam <- scanBam(system.file("extdata", "ex1.bam",
                               package="Rsamtools"))[[1]]

    checkIdentical(bam$flag, tbl[[2]])
    
    idx <- !is.na(bam$strand) & bam$strand=="+"
    checkIdentical(as.character(tbl[[10]][idx]),
                   as.character(bam$seq[idx]))
    checkIdentical(as.character(tbl[[11]][idx]),
                   as.character(bam$qual[idx]))

    idx <- !is.na(bam$strand) & bam$strand=="-"
    checkIdentical(as.character(tbl[[10]][idx]),
                   as.character(bam$seq[idx]))
    checkIdentical(as.character(tbl[[11]][idx]),
                   as.character(bam$qual[idx]))
    
}

test_scanBam_tag <- function()
{
    checkIdentical(character(0), bamTag(ScanBamParam()))
    tag <- c("MF", "Aq", "NM", "UQ", "H0", "H1")
    param <- ScanBamParam(tag=tag)
    checkTrue(validObject(param))
    checkIdentical(tag, bamTag(param))
    ## tags must be two letters
    checkException(ScanBamParam(tag="XYZ"), silent=TRUE)

    tag <- c("MF", "Aq", "NM", "UQ", "H0", "H1")
    param <- ScanBamParam(tag=tag)
    src <- system.file("unitTests", "cases", package="Rsamtools")
    fl <- file.path(src, "ex1_shuf1000.bam")
    bam <- scanBam(fl, param=param)[[1]][["tag"]]
    checkIdentical(tag, names(bam))
    checkTrue(all(1000L == sapply(bam, length)))
    exp <- structure(c(818L, 117L, 27L, 10L, 11L, 4L, 1L, 1L, 11L),
                     .Dim = 9L, .Dimnames = structure(list( c("0",
                     "1", "2", "3", "4", "5", "6", "7", NA)), .Names =
                     ""), class = "table")
    checkIdentical(exp, table(bam[["NM"]], useNA="always"))
}
