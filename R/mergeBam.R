setMethod(mergeBam, "character",
    function(files, destination, ..., region = GRanges(),
             overwrite = FALSE, header = character(), byQname = FALSE,
             addRG = FALSE, compressLevel1 = FALSE,
             indexDestination = FALSE)
{
    tryCatch({

        files <- sapply(files, .normalizePath)
        destination <- .normalizePath(destination)
        region <- local({
            x <- as(region, "GRanges")
            if (1L < length(x))
                stop("'region' must specify one range")
            sprintf("%s:%d-%d", as.character(seqnames(x)), start(x), end(x))
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
            indexBam(destination)

        destination

    }, error=function(err) {
        msg <- sprintf("'mergeBam' %s", conditionMessage(err))
        stop(msg)
    })
})
