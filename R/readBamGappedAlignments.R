setMethod(readBamGappedAlignments, "character",
          function(file, index=file, ..., which)
{
    if (missing(which))
        which <- RangesList()
    if (missing(index) && 0L == length(which))
        index <- character(0)
    bam <- open(BamFile(file, index), "rb")
    on.exit(close(bam))
    callGeneric(bam, character(), ..., which=which)
})
