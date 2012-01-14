setMethod(asBam, "character",
    function(file, destination, ...,
             overwrite=FALSE, indexDestination=TRUE)
{
    file <- .normalizePath(file)
    destination <- .normalizePath(destination)
    d0 <- paste(destination, "bam", sep=".")

    ofl <- tempfile()
    on.exit(unlink(ofl))
    tryCatch({
        if (!overwrite && file.exists(d0)) {
            msg <- sprintf("'%s' exists, '%s' is FALSE\n  %s: %s",
                           "destination", "overwrite", "destination",
                           d0)
            stop(msg)
        }
        result <- .Call(.as_bam, file, ofl)
        if (!file.exists(ofl))
            stop("failed to create 'BAM' file")
        if (indexDestination) {
            destination <- sortBam(ofl, destination)
            indexBam(destination)
        } else {
            destination <- d0
            .file.rename(ofl, destination)
        }
    }, error=function(err) {
        msg <- sprintf("'asBam' %s\n  SAM file: '%s'\n",
                       conditionMessage(err), file)
        stop(msg)
    })
    destination
})
