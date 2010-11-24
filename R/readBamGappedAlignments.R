setMethod(readBamGappedAlignments, "character",
          function(file, index, ..., which)
{
    if (missing(index))
        index <- .normalizePath(file)
    if (missing(which))
        which <- RangesList()
    bam <- open(BamFile(file, index), "rb")
    on.exit(close(bam))
    callGeneric(bam, character(), ..., which=which)
})
