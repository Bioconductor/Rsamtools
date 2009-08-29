fl <- system.file("extdata", "ex1.bam", package="Rsamtools")

test_countBam <- function()
{
    checkEquals(data.frame(space=NA, start=NA, end=NA, width=NA,
                           file=basename(fl), records=3307L,
                           nucleotides=116551L),
                countBam(fl))

    which <- RangesList(seq1=IRanges(1000, 2000),
                        seq2=IRanges(c(100, 1000), c(1000, 2000)))
    p1 <- ScanBamParam(which=which)
    exp <- cbind(as.data.frame(which),
                 file=basename(fl),
                 records=c(612L, 1168L, 642L),
                 nucleotides=c(21549L, 41200L, 22640L))
    checkIdentical(exp, countBam(fl, param=p1))

    which <- RangesList(seq2=IRanges(c(100, 1000), c(1000, 2000)),
                        seq1=IRanges(1000, 2000))
    p2 <- ScanBamParam(which=which)
    exp <- exp[c(2:3, 1),]
    rownames(exp) <- NULL
    checkIdentical(exp, countBam(fl, param=p2))
}
