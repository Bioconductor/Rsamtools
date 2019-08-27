test_io_check_exists <- function()
{
    .io_check_exists <- Rsamtools:::.io_check_exists

    fls <- c(
        tempfile(), tempfile(),
        "ftp://some.where/file",
        "http://some.where/file",
        "https://some.where/file",
        NA
    )
    file.create(fls[1])
    checkTrue(is.null(.io_check_exists(fls[1])))
    checkException(.io_check_exists(fls[2]), silent = TRUE)
    checkTrue(is.null(.io_check_exists(fls[3:6])))
}

test_catch_samtools <- function()
{
    fl <- system.file("unitTests", "cases", "ex1_unsort.bam",
                      package="Rsamtools")
    err <- warn <- FALSE
    tryCatch(suppressWarnings(withCallingHandlers({
        indexBam(fl)
    }, warning=function(msg) {
        warn <<- TRUE
    })), error=function(msg) {
        err <<- TRUE
    })
    checkTrue(isFALSE(warn))
    checkTrue(err)
}

test_catch_samtools_504 <- function()
{
    err <- FALSE
    tryCatch({
        scanBam("http://httpbin.org/status/504")
    }, error=function(err) {
        txt <- "failed to open BamFile:"
        err <<- startsWith(conditionMessage(err), txt)
    })
    checkTrue(err)
}

test_normalizePath <- function()
{
    .normalizePath <- Rsamtools:::.normalizePath
    fl <- tempfile()
    checkIdentical(fl, .normalizePath(fl))
    checkIdentical(fl, .normalizePath(factor(fl)))
}
