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
