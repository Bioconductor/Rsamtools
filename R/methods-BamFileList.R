BamFileList <-
    function(..., yieldSize=NA_integer_, obeyQname=FALSE)
{
    fls <- .RsamtoolsFileList(..., yieldSize=yieldSize, class="BamFile")
    if (!missing(obeyQname))
        obeyQname(fls) <- obeyQname
    fls
}

setMethod("summarizeOverlaps", c("GRanges", "BamFileList"),
    function(features, reads, mode, ignore.strand=FALSE, ..., 
             singleEnd=TRUE, param=ScanBamParam())
{
    .processBamFiles(features, reads, mode, ignore.strand, 
        ..., singleEnd=singleEnd, param=param)
})

setMethod("summarizeOverlaps", c("GRangesList", "BamFileList"),
    function(features, reads, mode, ignore.strand=FALSE, ..., 
             singleEnd=TRUE, param=ScanBamParam())
{
    .processBamFiles(features, reads, mode, ignore.strand, 
        ..., singleEnd=singleEnd, param=param)
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
