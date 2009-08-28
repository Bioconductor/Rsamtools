library(Rsamtools)
library(RUnit)

fl <- system.file("extdata", "ex1.bam", package="Rsamtools")

test_read_bam_header <- function()
{
    res <- .Call(Rsamtools:::.read_bam_header, fl, "rb", FALSE)
    checkIdentical(fl, names(res))
    exp <- structure(c(1575L, 1584L), .Names = c("seq1", "seq2"))
    checkIdentical(exp, res[[1]])
}

test_read_bam_header()
