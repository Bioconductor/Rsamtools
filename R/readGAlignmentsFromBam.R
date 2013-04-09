setMethod(readGAlignmentsFromBam, "character",
          function(file, index=file, ..., use.names=FALSE, param=NULL)
{
    if (missing(index) && (is.null(param) || 0L == length(bamWhich(param))))
        index <- character(0)
    bam <- open(BamFile(file, index), "rb")
    on.exit(close(bam))
    readGAlignmentsFromBam(bam, character(), ..., use.names=use.names,
                           param=param)
})

setMethod(readGappedReadsFromBam, "character",
          function(file, index=file, use.names=FALSE, param=NULL)
{
    if (missing(index) && (is.null(param) || 0L == length(bamWhich(param))))
        index <- character(0)
    bam <- open(BamFile(file, index), "rb")
    on.exit(close(bam))
    readGappedReadsFromBam(bam, character(), use.names=use.names, param=param)
})

setMethod(readGAlignmentPairsFromBam, "character",
          function(file, index=file, use.names=FALSE, param=NULL)
{
    if (missing(index) && (is.null(param) || 0L == length(bamWhich(param))))
        index <- character(0)
    bam <- open(BamFile(file, index), "rb")
    on.exit(close(bam))
    readGAlignmentPairsFromBam(bam, character(), use.names=use.names,
                               param=param)
})

setMethod(readGAlignmentsListFromBam, "character",
          function(file, index=file, ..., use.names=FALSE, 
                   param=ScanBamParam(), asProperPairs=TRUE)
{
    if (missing(index) && (is.null(param) || 0L == length(bamWhich(param))))
        index <- character(0)
    bam <- open(BamFile(file, index, obeyQname=TRUE), "rb")
    on.exit(close(bam))
    readGAlignmentsListFromBam(bam, character(), use.names=use.names,
                               param=param, asProperPairs=asProperPairs)
})

