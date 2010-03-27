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

setGeneric("filterBam",
           function(file, destination, index=file, ...)
           standardGeneric("filterBam"))

setGeneric("indexBam",
           function(files, ...) standardGeneric("indexBam"))

setGeneric("readPileup",
           function(file, ...) standardGeneric("readPileup"))

setGeneric("BamViews",
           function(bamPaths=character(0),
                    bamIndicies=bamPaths,
                    bamSamples=new("DataFrame", nrows=length(bamPaths),
                      rownames=make.unique(basename(bamPaths))),
                    bamRanges,
                    bamExperiment=list(), ...)
           standardGeneric("BamViews"),
           signature="bamRanges")

setGeneric("readBamGappedAlignments",
           function(file, index, ..., which)
           standardGeneric("readBamGappedAlignments"),
           signature="file")

