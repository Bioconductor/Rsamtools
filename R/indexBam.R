setMethod(indexBam, "character",
          function(files, ...)
{
    files <- .normalizePath(files)
    sapply(files, function(file) .Call(.index_bam, file))
})
