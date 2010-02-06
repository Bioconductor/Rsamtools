###

test_cigarToIRangesListByRName <- function()
{
    cigar <- c("30M5000N10M", "50M", "90M10M5I50M10D40M", "18M10I22M", "99I")
    rname <- c("chr6", "chr6", "chr2", "chr6", "chr2")
    pos <- c(101, 201, 1001,  301, 2001)
    ans <- as.list(cigarToIRangesListByRName(cigar=cigar, rname=rname,
                                             flag=NULL, pos=pos))
    ir2 <- IRanges(c(1001, 1091, 1101, 1161), c(1090, 1100, 1150, 1200))
    ir6 <- IRanges(c(101, 5131, 201, 301, 319), c(130, 5140, 250, 318, 340))
    ans0 <- list(chr2=ir2, chr6=ir6)
    checkEquals(ans, ans0)
}

