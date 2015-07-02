test_mapqFilter <- function() {
    fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
    p1 <- ScanBamParam(what="mapq")
    p2 <- ScanBamParam(what="mapq", mapqFilter=74)
    mapq1 <- scanBam(fl, param=p1)[[1]][["mapq"]]
    mapq2 <- scanBam(fl, param=p2)[[1]][["mapq"]]
    checkTrue(length(mapq1[!is.na(mapq1) & mapq1 == bamMapqFilter(p2)]) > 0)
    checkIdentical(mapq2, mapq1[!is.na(mapq1) & mapq1 >= bamMapqFilter(p2)])

    n <- countBam(fl, param=p2)
    checkIdentical(n$records, length(mapq2))

    checkException(ScanBamParam(mapqFilter=-1), silent=TRUE)
    checkException(ScanBamParam(mapqFilter=1:2), silent=TRUE)

    checkIdentical(bamMapqFilter(p2), 74L)
    bamMapqFilter(p2) <- 75
    checkIdentical(bamMapqFilter(p2), 75L)
    bamMapqFilter(p2) <- "76.1"
    checkIdentical(bamMapqFilter(p2), 76L)
}
