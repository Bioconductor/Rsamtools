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
           what="character",
           which="RangesList"))

setClass("Cigar",
         representation=representation("Sequence",
           .cigar="factor"),
         validity=.validity)
         
setClass("BamViews",
         representation=representation(
           bamPaths="character",
           bamIndicies="character",
           bamSamples="DataFrame",
           bamRanges="RangedData",
           bamExperiment="list"),
         validity=.validity)
           
### Since we are at a very early stage of prototyping this container, we'll
### call it Alignments0 for now.
### See http://wiki.fhcrc.org/bioc/Multiple_alignment_rep_v1 for the details
### of the class proposal.
setClass("Alignments0",
    representation(
        rname="factor",  # character factor
        strand="raw",
        cigar="character",
        ranges="CompressedNormalIRangesList"
    )
)
