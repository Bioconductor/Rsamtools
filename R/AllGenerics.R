setGeneric("ScanBamParam",
           function(flag=scanBamFlag(), simpleCigar=FALSE,
                    reverseComplement=FALSE, tag=character(0),
                    what=character(0), which)
           standardGeneric("ScanBamParam"),
           signature="which")

setGeneric("scanBam",
           function(file, index=file, ...,
                    param=ScanBamParam(what=scanBamWhat()))
           standardGeneric("scanBam"),
           signature="file")

setGeneric("scanBamHeader",
           function(files, ...) standardGeneric("scanBamHeader"))

setGeneric("countBam",
           function(file, index=file, ..., param=ScanBamParam())
           standardGeneric("countBam"),
           signature="file")
 
setGeneric("asBam",
           function(file, destination, ...)
           standardGeneric("asBam"))

setGeneric("sortBam",
           function(file, destination, ...)
           standardGeneric("sortBam"))

setGeneric("mergeBam",
           function(files, destination, ...)
           standardGeneric("mergeBam"),
           signature="files")

setGeneric("indexBam",
           function(files, ...) standardGeneric("indexBam"))

setGeneric("filterBam",
           function(file, destination, index=file, ...)
           standardGeneric("filterBam"))

setGeneric("quickCountBam",
           function(file, ..., param=ScanBamParam(), mainGroupsOnly=FALSE)
           standardGeneric("quickCountBam"),
           signature="file")

setGeneric("indexFa",
           function(file, ...) standardGeneric("indexFa"))

setGeneric("scanFaIndex",
           function(file, ...) standardGeneric("scanFaIndex"))

setGeneric("countFa",
           function(file, ...) standardGeneric("countFa"))

setGeneric("scanFa",
           function(file, param, ...) standardGeneric("scanFa"))

setGeneric("readPileup",
           function(file, ...) standardGeneric("readPileup"))

## bcf

setGeneric("ScanBcfParam",
           function(fixed=character(), info=character(), geno=character(), 
                    samples=character(), trimEmpty=TRUE, which, ...)
           standardGeneric("ScanBcfParam"),
           signature="which")

setGeneric("isOpen")

setGeneric("isIncomplete")

setGeneric("scanBcfHeader",
           function(file, ...) standardGeneric("scanBcfHeader"))

setGeneric("scanBcf",
           function(file, ...) standardGeneric("scanBcf"))

setGeneric("asBcf",
           function(file, dictionary, destination, ...,
                    overwrite=FALSE, indexDestination=TRUE)
           standardGeneric("asBcf"),
           signature="file")

setGeneric("indexBcf",
           function(file, ...) standardGeneric("indexBcf"))


## other

setGeneric("BamViews",
           function(bamPaths=character(0),
                    bamIndicies=bamPaths,
                    bamSamples=DataFrame(row.names=
                      make.unique(basename(bamPaths))),
                    bamRanges,
                    bamExperiment=list(), ...)
           standardGeneric("BamViews"),
           signature="bamRanges")

setGeneric("readGAlignmentsFromBam",
           function(file, index=file, ..., use.names=FALSE, param=NULL,
                    with.which_label=FALSE)
           standardGeneric("readGAlignmentsFromBam"),
           signature="file")

setGeneric("readGappedReadsFromBam",
           function(file, index=file, use.names=FALSE, param=NULL,
                    with.which_label=FALSE)
           standardGeneric("readGappedReadsFromBam"),
           signature="file")

setGeneric("readGAlignmentPairsFromBam",
           function(file, index=file, use.names=FALSE, param=NULL,
                    with.which_label=FALSE)
           standardGeneric("readGAlignmentPairsFromBam"),
           signature="file")

setGeneric("readGAlignmentsListFromBam",
           function(file, index=file, ..., use.names=FALSE, 
                    param=ScanBamParam(), with.which_label=FALSE)
           standardGeneric("readGAlignmentsListFromBam"),
           signature="file")

## tabix

setGeneric("seqnamesTabix", function(file, ...)
           standardGeneric("seqnamesTabix"))

setGeneric("headerTabix", function(file, ...)
           standardGeneric("headerTabix"))

setGeneric("scanTabix", function(file, ..., param)
           standardGeneric("scanTabix"))

setGeneric("yieldTabix", function(file, ..., yieldSize=1000000L)
           standardGeneric("yieldTabix"),
           signature="file")

## pileup

setGeneric(".asSpace", function(x) standardGeneric(".asSpace"))

setGeneric("applyPileups", function(files, FUN, ..., param)
           standardGeneric("applyPileups"),
           signature=c("files", "param"))

## RsamtoolsFile(s)

setGeneric("path",
           function(object, ...) standardGeneric("path"))

setGeneric("yieldSize",
           function(object, ...) standardGeneric("yieldSize"))

setGeneric("yieldSize<-",
           function(object, ..., value) standardGeneric("yieldSize<-"))

setGeneric("obeyQname",
           function(object, ...) standardGeneric("obeyQname"))

setGeneric("obeyQname<-",
           function(object, ..., value) standardGeneric("obeyQname<-"))

setGeneric("asMates",
           function(object, ...) standardGeneric("asMates"))

setGeneric("asMates<-",
           function(object, ..., value) standardGeneric("asMates<-"))

setGeneric("isOpen")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Old stuff.
###

readBamGappedAlignments <- function(...)
{
    .Deprecated("readGAlignmentsFromBam")
    readGAlignmentsFromBam(...)
}

readBamGappedReads <- function(...)
{
    .Deprecated("readGappedReadsFromBam")
    readGappedReadsFromBam(...)
}

readBamGappedAlignmentPairs <- function(...)
{
    .Deprecated("readGAlignmentPairsFromBam")
    readGAlignmentPairsFromBam(...)
}

readBamGAlignmentsList <- function(...)
{
    .Deprecated("readGAlignmentsListFromBam")
    readGAlignmentsListFromBam(...)
}

