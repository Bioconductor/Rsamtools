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
        txt <- "failed to open BamFile:"
        err <<- startsWith(conditionMessage(err), txt)
    })
    checkTrue(err)
}
