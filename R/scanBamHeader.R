setMethod(scanBamHeader, "character", 
          function(files, ...)
{
    files <- .normalizePath(files)
    lst <- lapply(files, function(file) {
        bam <- open(BamFile(file, character(0)))
        on.exit(close(bam))
        scanBamHeader(bam, ...)
    })
    names(lst) <- files
    lst
})
