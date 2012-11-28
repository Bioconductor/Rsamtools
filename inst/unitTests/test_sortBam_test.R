test_sortBam <- function() {
    fl0 <- system.file("extdata", "ex1.bam", package="Rsamtools")
    src <- system.file("unitTests", "cases", package="Rsamtools")
    fl <- file.path(src, "ex1_unsort.bam")
    ofl <- tempfile()
    checkTrue(file.create(ofl))
    on.exit(unlink(ofl))
    sorted <- sortBam(fl, ofl)
    exp <- scanBam(fl0)[[1]]
    obs <- scanBam(sorted)[[1]]
    checkIdentical(exp[["rname"]], obs[["rname"]])
    checkIdentical(Filter(Negate(is.na), exp[["pos"]]),
                   Filter(Negate(is.na), obs[["pos"]]))
}

test_sortBam_not_BAM_input <- function() {
    fl0 <- system.file("extdata", "ex1.sam", package="Rsamtools")
    checkException(sortBam(fl0, tempfile()), silent=TRUE)
}
