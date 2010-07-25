test_sortBam <- function() {
    fl0 <- system.file("extdata", "ex1.bam", package="Rsamtools")
    fl <- file.path("cases", "ex1_unsort.bam")
    sorted <- sortBam(fl, tempfile())
    exp <- scanBam(fl0)[[1]]
    obs <- scanBam(sorted)[[1]]
    checkIdentical(exp[["rname"]], obs[["rname"]])
    checkIdentical(Filter(Negate(is.na), exp[["pos"]]),
                   Filter(Negate(is.na), obs[["pos"]]))
}
