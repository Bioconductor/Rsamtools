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
      rname="factor", # character factor
      strand="raw",
      cigar="character",
      gapped_ranges="GappedRanges"
  )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some low-level helper functions for manipulating the @strand slot.
###

logicalAsCompactRawVector <- function(x)
{
    if (!is.logical(x))
        stop("'x' must be a logical vector")
    .Call(".logical_as_compact_raw_vector", x, PACKAGE="Rsamtools")
}

compactRawVectorAsLogical <- function(x, length.out)
{
    if (!is.raw(x))
        stop("'x' must be a raw vector")
    if (!isSingleNumber(length.out))
        stop("'length.out' must be a single number")
    if (!is.integer(length.out))
        length.out <- as.integer(length.out)
    .Call(".compact_raw_vector_as_logical", x, length.out, PACKAGE="Rsamtools")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor-like methods.
###

setMethod("length", "Alignments0", function(x) length(x@cigar))

setGeneric("rname", function(x) standardGeneric("rname"))
setMethod("rname", "Alignments0", function(x) x@rname)

setMethod("strand", "Alignments0",
    function(x)
    {
        is_minus <- compactRawVectorAsLogical(x@strand, length(x))
        ans <- factor(levels=c("+", "-", "*"))
        ans[!is_minus] <- "+"
        ans[is_minus] <- "-"
        return(ans)
    }
)

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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "coverage" method.
###

setMethod("coverage", "Alignments0",
    function (x, start=NA, end=NA, shift=0L, width=NULL, weight=1L, ...)
    {
        if (!identical(start, NA) || !identical(end, NA)
         || !identical(shift, 0L) || !is.null(width) || !identical(weight, 1L))
            stop("'start', 'end', 'shift', 'width' and 'weight' arguments ",
                 "are not supported yet, sorry!")
        irl <- cigarToIRangesList(cigar(x), rname(x), start(x))
        irl <- irl[elementLengths(irl) != 0]  # drop empty elements
        coverage(irl)
    }
)

