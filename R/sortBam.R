setMethod(sortBam, "character",
          function(file, destination, ...,
                   byQname=FALSE, maxMemory=512)
{
    file <- .normalizePath(file)
    destination <- .normalizePath(destination)
    result <- .Call(.sort_bam, file, destination, byQname,
                    as.integer(maxMemory))
    destination <- paste(result, "bam", sep=".")
    if (!file.exists(destination)) {
        msg <- sprintf("'sortBam' failed to create destination '%s'",
                       destination)
        stop(msg)
    }
    destination
})
