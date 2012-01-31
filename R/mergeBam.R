setMethod(mergeBam, "character",
    function(files, destination, ..., region = RangedData(),
             overwrite = FALSE, header = character(), byQname = FALSE,
             addRG = FALSE, compressLevel1 = FALSE,
             indexDestination = FALSE)
{
    tryCatch({

        files <- sapply(files, .normalizePath)
        destination <- .normalizePath(destination)
        region <- local({
            x <- as(region, "RangedData")
            if (1L < nrow(x))
                stop("'region' must specify one range")
            sprintf("%s:%d-%d", as.character(space(x)), start(x), end(x))
        })

        if (!overwrite && file.exists(destination)) {
            msg <- sprintf("'%s' exists, '%s' is FALSE\n  %s: %s",
                           "destination", "overwrite", "destination",
                           destination)
            stop(msg)
        }

        header <- .normalizePath(header)

        destination <-
            .Call(.merge_bam, files, destination, overwrite, header,
                  region, byQname, addRG, compressLevel1)
        if (indexDestination)
            destination <- indexBam(destination)

        destination

    }, error=function(err) {
        msg <- sprintf("'mergeBam' %s", conditionMessage(err))
        stop(msg)
    })
})
