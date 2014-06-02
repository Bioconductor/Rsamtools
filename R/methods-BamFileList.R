BamFileList <-
    function(..., yieldSize=NA_integer_, obeyQname=FALSE, asMates=FALSE,
             qnamePrefixEnd=NA, qnameSuffixStart=NA)
{
    fls <- .RsamtoolsFileList(..., yieldSize=yieldSize, class="BamFile")
    if (!missing(obeyQname))
        obeyQname(fls) <- obeyQname
    if (!missing(asMates))
        asMates(fls) <- asMates 
    if (!missing(qnamePrefixEnd))
        qnamePrefixEnd(fls) <- qnamePrefixEnd 
    if (!missing(qnameSuffixStart))
        qnameSuffixStart(fls) <- qnameSuffixStart
    fls
}

setMethod(obeyQname, "BamFileList",
    function(object, ...)
{
    sapply(object, obeyQname)
})

setReplaceMethod("obeyQname", "BamFileList", 
    function(object, ..., value)
{
    endoapply(object, `obeyQname<-`, value=value)
})

setMethod(asMates, "BamFileList",
    function(object, ...)
{
    sapply(object, asMates)
})

setReplaceMethod("asMates", "BamFileList", 
    function(object, ..., value)
{
    endoapply(object, `asMates<-`, value=value)
})

setMethod(qnamePrefixEnd, "BamFileList",
    function(object, ...)
{
    sapply(object, qnamePrefixEnd)
})

setReplaceMethod("qnamePrefixEnd", "BamFileList", 
    function(object, ..., value)
{
    endoapply(object, `qnamePrefixEnd<-`, value=value)
})

setMethod(qnameSuffixStart, "BamFileList",
    function(object, ...)
{
    sapply(object, qnameSuffixStart)
})

setReplaceMethod("qnameSuffixStart", "BamFileList", 
    function(object, ..., value)
{
    endoapply(object, `qnameSuffixStart<-`, value=value)
})

setMethod(seqinfo, "BamFileList",
    function(x)
{
    Reduce(merge, lapply(x, seqinfo))
})
