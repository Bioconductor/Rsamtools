BamFileList <-
    function(..., yieldSize=NA_integer_, obeyQname=FALSE, asMates=FALSE)
{
    fls <- .RsamtoolsFileList(..., yieldSize=yieldSize, class="BamFile")
    if (!missing(obeyQname))
        obeyQname(fls) <- obeyQname
    if (!missing(asMates))
        asMates(fls) <- asMates 
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
    endoapply(object, "obeyQname<-", value=value)
})

setMethod(asMates, "BamFileList",
    function(object, ...)
{
    sapply(object, asMates)
})

setReplaceMethod("asMates", "BamFileList", 
    function(object, ..., value)
{
    endoapply(object, "asMates<-", value=value)
})
