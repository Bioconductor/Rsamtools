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
           bamRanges="GRanges",
           bamExperiment="list"),
         validity=.validity)

setClass("ScanBcfParam",
    representation=representation(
      which="RangesList",
      info="character",
      geno="character",
      trimEmpty="logical"),
    prototype=prototype(
      trimEmpty=TRUE))

## RsamtoolsFile
.RsamtoolsFile_generator <- setRefClass("RsamtoolsFile",
    fields=list(.extptr="externalptr", path="character",
      index="character"))

.BamFile <- setRefClass("BamFile", contains="RsamtoolsFile")
.BcfFile <- setRefClass("BcfFile", contains="RsamtoolsFile",
   fields=list(mode="character"))
