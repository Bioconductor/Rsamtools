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
                    trimEmpty=TRUE, which, ...)
           standardGeneric("ScanBcfParam"),
           signature="which")

setGeneric("isOpen")

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

setGeneric("readBamGappedAlignments",
           function(file, index=file, use.names=FALSE, param=NULL)
           standardGeneric("readBamGappedAlignments"),
           signature="file")

setGeneric("readBamGappedReads",
           function(file, index=file, use.names=FALSE, param=NULL)
           standardGeneric("readBamGappedReads"),
           signature="file")

setGeneric("readBamGappedAlignmentPairs",
           function(file, index=file, use.names=FALSE, param=NULL)
           standardGeneric("readBamGappedAlignmentPairs"),
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

setGeneric("isOpen")
