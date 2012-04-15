setOldClass(c("bzfile", "connection"))
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

setClass("ScanBVcfParam",
    representation=representation(
      "VIRTUAL",
      which="RangesList",
      fixed="character",
      info="character",
      geno="character",
      trimEmpty="logical"),
    prototype=prototype(
      trimEmpty=TRUE))

setClass("ScanBcfParam", contains="ScanBVcfParam")

setClass("PileupParam",
    representation=representation(
      flag="integer",
      minBaseQuality="integer",
      minMapQuality="integer",
      minDepth="integer",
      maxDepth="integer",
      yieldSize="integer",
      yieldBy="character",
      yieldAll="logical",
      which="GRanges",
      what="character"),
    validity=.validity)

## RsamtoolsFile(s)
.RsamtoolsFile_generator <- setRefClass("RsamtoolsFile",
    fields=list(.extptr="externalptr", path="character",
      index="character", yieldSize="integer"))

.BamFile <- setRefClass("BamFile", contains="RsamtoolsFile")
.BcfFile <- setRefClass("BcfFile", contains="RsamtoolsFile",
   fields=list(mode="character"))
.TabixFile <- setRefClass("TabixFile", contains="RsamtoolsFile")
.FaFile <- setRefClass("FaFile", contains="RsamtoolsFile")

setClass("RsamtoolsFileList", contains=c("SimpleList", "VIRTUAL"))

setClass("BamFileList", contains="RsamtoolsFileList",
         prototype=prototype(elementType="BamFile"))

setClass("BcfFileList", contains="RsamtoolsFileList",
         prototype=prototype(elementType="BcfFile"))

setClass("TabixFileList", contains="RsamtoolsFileList",
         prototype=prototype(elementType="TabixFile"))

setClass("FaFileList", contains="RsamtoolsFileList",
         prototype=prototype(elementType="FaFile"))

setClass("PileupFiles", contains="BamFileList",
         representation=representation(param="PileupParam"))
