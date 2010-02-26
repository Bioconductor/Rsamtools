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

### Incomplete implementation (the genomic ranges of the alignments are
### missing). The full implementation is achieved by each concrete class
### below. See inst/doc/GappedAlignments-timings.txt for the timings (and
### memory footprint) of the various GappedAlignments implementations.
setClass("GappedAlignments",
    contains="Sequence",
    representation(
        "VIRTUAL",
        cigar="character"             # extended CIGAR (see SAM format specs)
        #mismatches="characterORNULL", # see MD optional field in SAM format specs
        #values="DataFrame"
    )
)

### First GappedAlignments implementation: Alignments0
### The genomic ranges of the alignments are stored in 3 slots: 'rname',
### 'strand' and 'rglist'.
### See http://wiki.fhcrc.org/bioc/Multiple_alignment_rep_v1 for the details
### of the class proposal.
setClass("Alignments0",
    contains="GappedAlignments",
    representation(
        rname="factor",               # character factor
        strand="raw",
        rglist="CompressedNormalIRangesList"
    )
)

### Second GappedAlignments implementation: Alignments1
### The genomic ranges of the alignments are stored in the 'grglist' slot of
### type GRangesList.
setClass("Alignments1",
    contains="GappedAlignments",
    representation(
        grglist="GRangesList"
    )
)

### Third GappedAlignments implementation: Alignments2
### The genomic ranges of the alignments are stored in 3 slots: 'rname',
### 'strand' and 'start'.
setClass("Alignments2",
    contains="GappedAlignments",
    representation(
        rname="factor",               # character factor
        strand="raw",
        start="integer"               # POS field in SAM
    )
)

