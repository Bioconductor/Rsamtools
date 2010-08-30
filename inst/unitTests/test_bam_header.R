fl <- system.file("extdata", "ex1.bam", package="Rsamtools")

test_readBamHeader <- function()
{
    res <- scanBamHeader(fl)
    checkIdentical(Rsamtools:::.normalizePath(fl), names(res))
    exp <- structure(c(1575L, 1584L), .Names = c("seq1", "seq2"))
    checkIdentical(exp, res[[1]][["targets"]])
}
