setGeneric("ScanBamParam",
           function(flag=scanBamFlag(), simpleCigar=FALSE,
                    reverseComplement=FALSE, tag=character(0),
                    what=scanBamWhat(), which)
           standardGeneric("ScanBamParam"),
           signature="which")

setGeneric("scanBam",
           function(file, index=file, ..., param=ScanBamParam())
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
           function(file, param=GRanges(), ...)
           standardGeneric("scanFa"))

setGeneric("readPileup",
           function(file, ...) standardGeneric("readPileup"))

## bcf

setGeneric("ScanBcfParam",
           function(info=character(), geno=character(), trimEmpty=TRUE,
                    which)
           standardGeneric("ScanBcfParam"),
           signature="which")

setGeneric("isOpen")

setGeneric("scanBcfHeader",
           function(file, ...) standardGeneric("scanBcfHeader"))

setGeneric("scanBcf",
           function(file, ...) standardGeneric("scanBcf"))

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
           function(file, index, ..., which)
           standardGeneric("readBamGappedAlignments"),
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

## RsamtoolsFile

setGeneric("isOpen")
