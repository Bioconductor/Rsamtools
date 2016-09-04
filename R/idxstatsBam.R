setMethod(idxstatsBam, "character",
    function(file, index=file, ...)
{
    index <- .normalizePath(index)
    bam <- open(BamFile(file, index), "rb")
    on.exit(close(bam))
    idxstatsBam(bam, ...)
})
