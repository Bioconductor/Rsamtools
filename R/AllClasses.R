setOldClass("connection")
#setOldClass(c("file", "connection"))
#setOldClass(c("url", "connection"))
#setOldClass(c("gzfile", "connection"))
#setOldClass(c("bzfile", "connection"))
#setOldClass(c("unz", "connection"))
#setOldClass(c("file", "connection"))
#setOldClass(c("pipe", "connection"))
#setOldClass(c("fifo", "connection"))

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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### GappedAlignments objects
###

setClass("GappedAlignments", contains="Sequence", representation("VIRTUAL"))

### First GappedAlignments implementation: Alignments0
### See http://wiki.fhcrc.org/bioc/Multiple_alignment_rep_v1 for the details
### of the class proposal.
setClass("Alignments0",
    contains="GappedAlignments",
    representation(
        rname="factor",               # character factor
        strand="raw",
        cigar="character",            # extended CIGAR (see SAM format specs)
        ranges="CompressedNormalIRangesList"
    )
)

### Second GappedAlignments implementation: Alignments1
### The implementation of Alignments1 is based on the new GenomicFeatureList
### container defined in BSgenome. With just the 'cigar' slot as a proper
### additional slot, Alignments1 is equivalent to Alignments0 i.e. it allows
### storing the same information than Alignments0 (with some internal
### redundancy though that could make it slightly bigger).
setClass("Alignments1",
    contains=c("GappedAlignments", "GenomicFeatureList"),
    representation(
        cigar="character"             # extended CIGAR (see SAM format specs)
        #mismatches="characterORNULL", # see MD optional field in SAM format specs
        #values="DataFrame"
    ),
    prototype(elementType="GenomicFeature")
)

