test_asBam <- function()
{
    src <- system.file("unitTests", "cases", package="Rsamtools")
    fl <- file.path(src, "ex1.sam.gz")
    ofl <- tempfile()
    on.exit(unlink(ofl))
    bam <- asBam(fl, ofl)
    checkIdentical(bam, paste(ofl, "bam", sep="."))

    which <- GRanges("seq2", IRanges(1000, 2000))
    res <- scanBam(bam, param=ScanBamParam(which=which, what="rname"))[[1]]
    checkIdentical(642L, length(res[["rname"]]))
    checkIdentical("seq2", as.character(unique(res[["rname"]])))

    checkException(asBam(fl, ofl), silent=TRUE)
}
