setOldClass("connection")
setOldClass(c("file", "connection"))
setOldClass(c("url", "connection"))
setOldClass(c("gzfile", "connection"))
setOldClass(c("bzfile", "connection"))
setOldClass(c("unz", "connection"))
setOldClass(c("file", "connection"))
setOldClass(c("pipe", "connection"))
setOldClass(c("fifo", "connection"))

setGeneric(".validity", function(object) standardGeneric(".validity"))

setClass("ScanBamParam",
         representation=representation(
           flag="integer",
           simpleCigar="logical",
           reverseComplement="logical",
           tag="character",
           what="character",
           which="RangesList"))
         
setClass("BamViews",
         representation=representation(
           bamPaths="character",
           bamIndicies="character",
           bamSamples="DataFrame",
           bamRanges="RangedData",
           bamExperiment="list"),
         validity=.validity)

