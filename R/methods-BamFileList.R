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

setMethod("summarizeOverlaps", c("GRanges", "BamFileList"),
    function(features, reads, mode, ignore.strand=FALSE, ..., 
             inter.feature=TRUE, singleEnd=TRUE, fragments=TRUE,
             param=ScanBamParam())
{
    .dispatchBamFiles(features, reads, mode, ignore.strand, ..., 
                      inter.feature=inter.feature, singleEnd=singleEnd, 
                      fragments=fragments, param=param)
})

setMethod("summarizeOverlaps", c("GRangesList", "BamFileList"),
    function(features, reads, mode, ignore.strand=FALSE, ..., 
             inter.feature=TRUE, singleEnd=TRUE, fragments=TRUE,
             param=ScanBamParam())
{
    .dispatchBamFiles(features, reads, mode, ignore.strand, ..., 
                      inter.feature=inter.feature, singleEnd=singleEnd, 
                      fragments=fragments, param=param)
})

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
