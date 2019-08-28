setMethod(idxstatsBam, "character",
    function(file, index, ...)
{
    index <- if (missing(index))
                 character(0)
             else .normalizePath(index)
    bam <- open(BamFile(file, index), "rb")
    on.exit(close(bam))
    idxstatsBam(bam, ...)
})
