library(Rsamtools)
library(RUnit)

fl <- system.file("extdata", "ex1.bam", package="Rsamtools")

test_countBam <- function()
{
    checkEquals(3307L, countBam(fl))

    which <- RangesList(seq1=IRanges(1000, 2000),
                        seq2=IRanges(c(100, 1000), c(1000, 2000)))
    p1 <- ScanBamParam(which=which)
    vals <- c(612L, 1168L, 642L)
    names(vals) <- 
        paste(rep(names(which), sapply(which, length)), ":",
              start(which), "-", end(which), sep="")
    checkIdentical(vals, countBam(fl, param=p1))
}





