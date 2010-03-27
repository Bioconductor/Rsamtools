setMethod(readBamGappedAlignments, "character",
          function(file, index, ..., which)
{
    ans.subtype="Alignments2"
    if (missing(index))
        index <- file
    if (missing(which))
        which <- RangesList()
    param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE,
                                           isDuplicate=FALSE),
                          what=c("rname", "strand", "pos", "cigar"),
                          which=which)
    bam <- scanBam(file, index=index, param=param)
    ## unlist(list(factor())) returns integer(0), so exit early if all
    ## values are empty
    bam <- bam[sapply(bam, function(x) length(x$rname) != 0)]
    if (0L == length(bam))
        return(get(ans.subtype)())
    rname <- unlist(unname(lapply(bam, "[[", "rname")))
    strand <- unlist(unname(lapply(bam, "[[", "strand")))
    pos <- unlist(unname(lapply(bam, "[[", "pos")))
    cigar <- unlist(unname(lapply(bam, "[[", "cigar")))
    ## Calls the appropriate constructor.
    get(ans.subtype)(rname=rname, strand=strand, pos=pos, cigar=cigar)
})
