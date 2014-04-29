fl <- system.file("extdata", "ex1.bam", package="Rsamtools")

as_data.frame <- function(x) {
    data.frame(space=factor(names(unlist(x)),
                 levels=seqlevels(x)),
               start=start(unlist(x)),
               end=end(unlist(x)),
               width=width(unlist(x)))
}

test_countBam <- function()
{
    checkEquals(data.frame(space=NA, start=NA, end=NA, width=NA,
                           file=basename(fl), records=3307L,
                           nucleotides=116551L),
                countBam(fl))

    which <- RangesList(seq1=IRanges(1000, 2000),
                        seq2=IRanges(c(100, 1000), c(1000, 2000)))
    p1 <- ScanBamParam(which=which)
    exp <- cbind(as_data.frame(which),
                 file=basename(fl),
                 records=c(612L, 1169L, 642L),
                 nucleotides=c(21549, 41235, 22640))
    checkIdentical(exp, countBam(fl, param=p1))

    which <- RangesList(seq2=IRanges(c(100, 1000), c(1000, 2000)),
                        seq1=IRanges(1000, 2000))
    p2 <- ScanBamParam(which=which)
    exp <- merge(as_data.frame(which), exp, sort=FALSE)
    rownames(exp) <- NULL
    checkIdentical(exp, countBam(fl, param=p2))
}

test_countBam_index <- function()
{
    which <- RangesList(seq1=IRanges(1000, 2000),
                        seq2=IRanges(c(100, 1000), c(1000, 2000)))
    p1 <- ScanBamParam(which=which)
    exp <- cbind(as_data.frame(which),
                 file="ex1_noindex.bam",
                 records=c(612L, 1169L, 642L),
                 nucleotides=c(21549, 41235, 22640))

    src <- system.file("unitTests", "cases", package="Rsamtools")
    fl <- file.path(src, "ex1_noindex.bam")
    idx <- system.file("extdata", "ex1.bam", package="Rsamtools")
    checkIdentical(exp, countBam(fl, idx, param=p1))

    checkException({
        suppressWarnings(countBam(fl, tempfile(), param=p1))
    }, silent=TRUE)
}
