setMethod(asSam, "character",
    function(file, destination=sub("\\.bam", "", file), ..., overwrite=FALSE)
{
    file <- .normalizePath(file)
    destination <- .normalizePath(destination)
    d0 <- paste(destination, "sam", sep=".")

    if (!overwrite && file.exists(d0)) {
        msg <- sprintf("'%s' exists, '%s' is FALSE\n  %s: %s",
                       "destination", "overwrite", "destination",
                       d0)
        stop(msg)
    }

    tryCatch({
        result <- .Call(.as_bam, file, d0, FALSE)
        if (!file.exists(d0))
            stop("failed to create 'SAM' file")
    }, error=function(err) {
        msg <- sprintf("'asSam' %s\n  SAM file: '%s'\n",
                       conditionMessage(err), file)
        stop(msg)
    })
    d0
})

setMethod(asBam, "character",
    function(file, destination=sub("\\.sam(\\.gz)?", "", file), ...,
             overwrite=FALSE, indexDestination=TRUE)
{
    file <- .normalizePath(file)
    destination <- .normalizePath(destination)
    d0 <- paste(destination, "bam", sep=".")

    ofl <- tempfile()
    on.exit(unlink(ofl))
    if (!overwrite && file.exists(d0)) {
        msg <- sprintf("'%s' exists, '%s' is FALSE\n  %s: %s",
                       "destination", "overwrite", "destination",
                       d0)
        stop(msg)
    }
    tryCatch({
        result <- .Call(.as_bam, file, ofl, TRUE)
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
