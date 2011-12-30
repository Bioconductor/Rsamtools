test_unpackVcfInfo <- function()
{
    .unpackVcfInfo <- Rsamtools:::.unpackVcfInfo
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")

    ## info
    fmt <- scanVcfHeader(fl)[[1]][["Header"]][["INFO"]]
    info <- scanVcf(fl)[[1]]$INFO
    res <- .unpackVcfInfo(info, rownames(fmt), fmt$Number, fmt$Type)

    checkIdentical(as.integer(c(3, 3, 2, 3, 3)), res$NS)
    checkIdentical(as.integer(c(14, 11, 10, 13, 9)), res$DP)
    checkEquals(list(.5, .017, c(.333, .667), NA_real_, NA_real_),
                res$AF)
    checkIdentical(c(TRUE, FALSE, FALSE, FALSE, FALSE), res$DB)
    checkIdentical(rep(FALSE, 5), res$H2)
}

test_unpackVcfGeno <- function()
{
    .unpackVcfGeno <- Rsamtools:::.unpackVcfGeno
    fl <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")

    ## geno
    fmt <- scanVcfHeader(fl)[[1]][["Header"]][["FORMAT"]]
    geno <- scanVcf(fl)[[1]]$GENO
    res <- .unpackVcfGeno(geno, rownames(fmt), fmt$Number, fmt$Type)

    checkEquals(typeof(unlist(res$GT)), "character")
    checkIdentical(lapply(res, class), list(GT="matrix", GQ="matrix",
                   DP="matrix", HQ="array"))
    checkIdentical(as.integer(c(1, 3, 6, 7, 4, 8, 5, 0, 4, 2)),
                   as.vector(res$DP))
} 


