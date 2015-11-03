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
    checkTrue(warn)
    checkTrue(err)
}

test_catch_samtools_504 <- function()
{
    err <- FALSE
    tryCatch({
        scanBam("http://httpbin.org/status/504")
    }, error=function(err) {
        txt <- "failed to open BamFile: [khttp_connect_file] fail to open file (HTTP code: 504).\n"
        err <<- identical(conditionMessage(err), txt)
    })
    checkTrue(err)
}
