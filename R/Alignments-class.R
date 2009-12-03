### =========================================================================
### Alignments objects
### -------------------------------------------------------------------------
###

### Since we are at a very early stage of prototyping this container, we'll
### call it Alignments0 for now.
### See http://wiki.fhcrc.org/bioc/Multiple_alignment_rep_v1 for the details
### of the class proposal.
setClass("Alignments0",
  contains="Ranges",
  representation(
      rname="factor",     # character factor
      strand="raw",
      cigar="character",
      gapped_ranges="GappedRanges"
  )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some helper function for @strand slot manipulation.
###

logicalAsCompactRawVector <- function(x)
    .Call(".logical_as_compact_raw_vector", x, PACKAGE="Rsamtools")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor-like methods.
###

setMethod("length", "Alignments0", function(x) length(x@cigar))

setGeneric("rname", function(x) standardGeneric("rname"))
setMethod("rname", "Alignments0", function(x) x@rname)

### FIXME: Return a character factor with levels + - *
setMethod("strand", "Alignments0", function(x) x@strand)

setGeneric("cigar", function(x) standardGeneric("cigar"))
setMethod("cigar", "Alignments0", function(x) x@cigar)

setGeneric("gappedRanges", function(x) standardGeneric("gappedRanges"))
setMethod("gappedRanges", "Alignments0", function(x) x@gapped_ranges)

setMethod("start", "Alignments0", function(x, ...) start(gappedRanges(x)))

setMethod("end", "Alignments0", function(x, ...) end(gappedRanges(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

### IMPLEMENT ME!


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

### IMPLEMENT ME! (we just inherit the method for Ranges object right now)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors.
###

### This is our only constructor for now.
readBAMasAligments <- function(file, index=file, which=RangesList())
{
    param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE,
                                           isDuplicate=FALSE),
                          what=c("rname", "strand", "pos", "cigar"),
                          which=which)
    bam <- scanBam(file, index=index, param=param)[[1]]
    ans_rname <- bam$rname
    if (!is.factor(ans_rname))
        ans_rname <- as.factor(ans_rname)
    ans_strand <- logicalAsCompactRawVector(bam$strand == "-")
    ans_cigar <- bam$cigar@.cigar
    if (is.factor(ans_cigar))
        ans_cigar <- as.vector(ans_cigar)
    ans_gapped_ranges <- cigarToGappedRanges(ans_cigar, bam$pos)
    new("Alignments0", rname=ans_rname, strand=ans_strand,
                       cigar=ans_cigar, gapped_ranges=ans_gapped_ranges)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setAs("Alignments0", "GappedRanges", function(from) gappedRanges(from))
setAs("Alignments0", "CompressedIRangesList",
    function(from) as(gappedRanges(from), "CompressedIRangesList")
)
setAs("Alignments0", "IRangesList",
    function(from) as(gappedRanges(from), "IRangesList")
)
setAs("Alignments0", "RangesList",
    function(from) as(gappedRanges(from), "RangesList")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###


