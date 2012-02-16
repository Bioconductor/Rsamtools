setMethod(readBamGappedAlignments, "character",
          function(file, index=file, use.names=FALSE, param=NULL)
{
    if (missing(index) && (is.null(param) || 0L == length(bamWhich(param))))
        index <- character(0)
    bam <- open(BamFile(file, index), "rb")
    on.exit(close(bam))
    callGeneric(bam, character(), use.names=use.names, param=param)
})

setMethod(readBamGappedReads, "character",
          function(file, index=file, use.names=FALSE, param=NULL)
{
    if (missing(index) && (is.null(param) || 0L == length(bamWhich(param))))
        index <- character(0)
    bam <- open(BamFile(file, index), "rb")
    on.exit(close(bam))
    callGeneric(bam, character(), use.names=use.names, param=param)
})

setMethod(readBamGappedAlignmentPairs, "character",
          function(file, index=file, use.names=FALSE, param=NULL)
{
    if (missing(index) && (is.null(param) || 0L == length(bamWhich(param))))
        index <- character(0)
    bam <- open(BamFile(file, index), "rb")
    on.exit(close(bam))
    callGeneric(bam, character(), use.names=use.names, param=param)
})

