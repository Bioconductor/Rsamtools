setMethod(readBamGappedAlignments, "character",
          function(file, index, ..., which)
{
    if (missing(index)) {
        tmp <- paste(.normalizePath(file), "bai", sep=".")
        index <- if (file.exists(tmp)) file else character(0)
    }
    if (missing(which))
        which <- RangesList()
    bam <- open(BamFile(file, index), "rb")
    on.exit(close(bam))
    callGeneric(bam, character(), ..., which=which)
})
