setGeneric("ScanBamParam",
           function(flag=scanBamFlag(), simpleCigar=FALSE,
                    reverseComplement=FALSE, tag=character(0),
                    tagFilter=list(), what=character(0), which,
                    mapqFilter=NA_integer_)
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
 
setGeneric("idxstatsBam",
           function(file, index=file, ...) standardGeneric("idxstatsBam"),
           signature="file")

setGeneric("asBam",
           function(file, destination = sub("\\.sam(\\.gz)?", "", file), ...)
           standardGeneric("asBam"))

setGeneric("asSam",
           function(file, destination = sub("\\.bam", "", file), ...)
           standardGeneric("asSam"))

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

setGeneric("quickBamFlagSummary",
           function(file, ..., param=ScanBamParam(), main.groups.only=FALSE)
           standardGeneric("quickBamFlagSummary"),
           signature="file")

setGeneric("indexFa",
           function(file, ...) standardGeneric("indexFa"))

setGeneric("scanFaIndex",
           function(file, ...) standardGeneric("scanFaIndex"))

setGeneric("countFa",
           function(file, ...) standardGeneric("countFa"))

setGeneric("scanFa",
           function(file, param, ...,
                    as=c("DNAStringSet", "RNAStringSet", "AAStringSet"))
           standardGeneric("scanFa"),
           signature=c("file", "param"))

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

setGeneric("pileup",
    function(file, index=file, ..., scanBamParam=ScanBamParam(),
             pileupParam=PileupParam())
        standardGeneric("pileup"),
    signature=signature("file"))

setGeneric(".asRegions", function(x) standardGeneric(".asRegions"))

setGeneric("applyPileups",
           function(files, FUN, ..., param)
           {
               .Deprecated("pileup")
               standardGeneric("applyPileups")
           },
           signature=c("files", "param"))

## RsamtoolsFile(s)

setGeneric("index",
           function(object, ..., asNA=TRUE)
           standardGeneric("index"),
           signature="object")

setGeneric("index<-",
           function(object, ..., value) standardGeneric("index<-"))

setGeneric("gzindex",
           function(object, asNA=TRUE)
           standardGeneric("gzindex"),
           signature="object")

setGeneric("gzindex<-",
           function(object, value) standardGeneric("gzindex<-"))

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

setGeneric("qnamePrefixEnd",
           function(object, ...) standardGeneric("qnamePrefixEnd"))

setGeneric("qnamePrefixEnd<-",
           function(object, ..., value) standardGeneric("qnamePrefixEnd<-"))

setGeneric("qnameSuffixStart",
           function(object, ...) standardGeneric("qnameSuffixStart"))

setGeneric("qnameSuffixStart<-",
           function(object, ..., value) standardGeneric("qnameSuffixStart<-"))

setGeneric("isOpen")

setGeneric("testPairedEndBam", function(file, index=file, ...) 
           standardGeneric("testPairedEndBam"),
           signature="file")
