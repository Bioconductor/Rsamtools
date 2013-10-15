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

.summarizeOverlaps_character <-
    function(features, reads, mode, ignore.strand=FALSE, ...,
        yieldSize=1000000L, inter.feature=TRUE, singleEnd=TRUE,
        fragments=TRUE, param=ScanBamParam())
{
    if (is.null(names(reads))) {
        if (any(duplicated(reads)))
            stop("duplicate 'reads' paths not allowed; use distinct names()")
    } else if (any(duplicated(names(reads))))
        stop("duplicate 'names(reads)' file paths not allowed")
    reads <- BamFileList(reads, yieldSize=yieldSize, obeyQname=FALSE,
                         asMates=!singleEnd)
    summarizeOverlaps(features, reads, mode, ignore.strand, ...,
                      inter.feature=inter.feature, singleEnd=singleEnd,
                      fragments=fragments)
}

setMethod("summarizeOverlaps", c("GRanges", "character"),
    .summarizeOverlaps_character)

setMethod("summarizeOverlaps", c("GRangesList", "character"),
    .summarizeOverlaps_character)

.summarizeOverlaps_BamFileList <-
    function(features, reads, mode, ignore.strand=FALSE, ..., 
             inter.feature=TRUE, singleEnd=TRUE, fragments=TRUE,
             param=ScanBamParam())
{
    if (any(duplicated(names(reads))))
        stop("duplicate 'names(reads)' not allowed")
    .checkArgs(reads, singleEnd, fragments)
    .dispatchBamFiles(features, reads, mode, ignore.strand, ..., 
                      inter.feature=inter.feature, singleEnd=singleEnd, 
                      fragments=fragments, param=param)
}

setMethod("summarizeOverlaps", c("GRanges", "BamFileList"),
    .summarizeOverlaps_BamFileList)

setMethod("summarizeOverlaps", c("GRangesList", "BamFileList"),
    .summarizeOverlaps_BamFileList)

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
