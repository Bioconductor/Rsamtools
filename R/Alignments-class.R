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
        rname="factor",  # character factor
        strand="raw",
        cigar="character",
        gapped_ranges="GappedRanges"
    ),
    prototype(
        elementType="NormalIRanges"
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

subsetCompactRawVector <- function(x, i)
{
    if (!is.raw(x))
        stop("'x' must be a raw vector")
    if (!is.integer(i))
        stop("'i' must be an integer vector")
    .Call(".subset_compact_raw_vector", x, i, PACKAGE="Rsamtools")
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

setGeneric("qwidth", function(x) standardGeneric("qwidth"))
setMethod("qwidth", "Alignments0", function(x) cigarToQWidth(cigar(x)))

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

setMethod("as.data.frame", "Alignments0",
    function(x, row.names=NULL, optional=FALSE, ...)
    {
        if (!(is.null(row.names) || is.character(row.names)))
            stop("'row.names' must be NULL or a character vector")
        ans <- data.frame(rname=rname(x),
                          strand=strand(x),
                          cigar=cigar(x),
                          start=start(x),
                          end=end(x),
                          width=width(x),
                          row.names=row.names,
                          check.rows=TRUE,
                          check.names=FALSE,
                          stringsAsFactors=FALSE)
        return(ans)
    }
)

setMethod("show", "Alignments0",
    function(object)
    {
        lo <- length(object)
        cat(class(object), " of length ", lo, "\n", sep="")
        if (lo == 0L) {
            return(NULL)
        } else if (lo < 20L) {
            showme <-
              as.data.frame(object,
                            row.names=paste("[", seq_len(lo), "]", sep=""))
        } else {
            ## Use of as.vector() here is to prevent c() to do silly things
            ## when 'x' is a factor! (Try 'c(factor(LETTERS))', yes it's
            ## documented that c() will drop the attributes but still, this
            ## doesn't make sense).
            sketch <- function(x)
                          c(as.vector(x[1:9]),
                            "...",
                            as.vector(x[(length(x)-8L):length(x)]))
            showme <-
              data.frame(rname=sketch(rname(object)),
                         strand=sketch(strand(object)),
                         cigar=sketch(cigar(object)),
                         start=sketch(start(object)),
                         end=sketch(end(object)),
                         width=sketch(width(object)),
                         row.names=c(paste("[", 1:9, "]", sep=""), "...",
                                     paste("[", (lo-8L):lo, "]", sep="")),
                         check.rows=TRUE,
                         check.names=FALSE,
                         stringsAsFactors=FALSE)
        }
        show(showme)
    }
)


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

### Supported 'i' types: numeric vector, logical vector, NULL and missing.
setMethod("[", "Alignments0",
    function(x, i, j, ... , drop=TRUE)
    {
        if (!missing(j) || length(list(...)) > 0L)
            stop("invalid subsetting")
        if (missing(i))
            return(x)
        if (!is.atomic(i))
            stop("invalid subscript type")
        lx <- length(x)
        if (is.numeric(i)) {
            if (min(i) < 0L)
                i <- seq_len(lx)[i]
            else if (!is.integer(i))
                i <- as.integer(i)
        } else if (is.logical(i)) {
            if (length(i) > lx)
                stop("subscript out of bounds")
            i <- seq_len(lx)[i]
        } else if (!is.null(i)) {
            stop("invalid subscript type")
        }
        x@strand <- subsetCompactRawVector(x@strand, i)
        x@rname <- x@rname[i]
        x@cigar <- x@cigar[i]
        x@gapped_ranges <- x@gapped_ranges[i]
        x
    }
)

setMethod("[[", "Alignments0",
    function(x, i, j, ..., exact=TRUE) gappedRanges(x)[[i]]
)

### Without this definition, we inherit the method for Sequence objects
### which returns the same thing but is thousands of times slower!
setMethod("elementLengths", "Alignments0",
    function(x) elementLengths(gappedRanges(x))
)


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

