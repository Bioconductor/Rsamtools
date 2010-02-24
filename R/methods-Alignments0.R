### =========================================================================
### Alignments0 objects
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some low-level helper functions for manipulating the @strand slot.
###

.logicalAsCompactRawVector <- function(x)
{
    if (!is.logical(x))
        stop("'x' must be a logical vector")
    .Call(".logical_as_compact_raw_vector", x, PACKAGE="Rsamtools")
}

.compactRawVectorAsLogical <- function(x, length.out)
{
    if (!is.raw(x))
        stop("'x' must be a raw vector")
    if (!isSingleNumber(length.out))
        stop("'length.out' must be a single number")
    if (!is.integer(length.out))
        length.out <- as.integer(length.out)
    .Call(".compact_raw_vector_as_logical", x, length.out, PACKAGE="Rsamtools")
}

.subsetCompactRawVector <- function(x, i)
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

setMethod("rname", "Alignments0", function(x) x@rname)

setMethod("strand", "Alignments0",
    function(x) strand(.compactRawVectorAsLogical(x@strand, length(x)))
)

setMethod("granges", "Alignments0",
    function(x)
        GappedAlignmentsAsGRangesList(rname(x), strand(x), ranges(x))
)

setMethod("ranges", "Alignments0", function(x) x@ranges)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors.
###

Alignments0 <- function(rname=factor(), strand=BSgenome::strand(),
                        pos=integer(0), cigar=character(0))
{
    if (!is.factor(rname) || !is.character(levels(rname))) {
        if (!is.character(rname))
            stop("'rname' must be a character vector/factor")
        rname <- as.factor(rname)
    }
    if (any(is.na(rname)))
        stop("'rname' cannot have NAs")
    if (!is.factor(strand) || !identical(levels(strand), .STRAND_LEVELS))
        stop("'strand' must be a character factor")
    if (!is.character(cigar) || any(is.na(cigar)))
        stop("'cigar' must be a character vector with no NAs")
    ranges <- cigarToIRangesListByAlignment(cigar, pos)
    strand <- .logicalAsCompactRawVector(strand == "-")
    new("Alignments0", rname=rname, strand=strand,
                       cigar=cigar, ranges=ranges)
}

### This is our only constructor for now.
setMethod(readBAMasAlignments0, "character", 
          function(file, index, ..., which)
{
    if (missing(index))
        index <- file
    if (missing(which))
        which <- RangesList()
    param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE,
                                           isDuplicate=FALSE),
                          what=c("rname", "strand", "pos", "cigar"),
                          which=which)
    bam <- scanBam(file, index=index, param=param)
    ## unlist(list(factor())) returns integer(0), so exit early if all
    ## values are empty
    if (all(sapply(bam, function(x) length(x$rname) == 0)))
        return(Alignments0())
    rname <- unlist(unname(lapply(bam, "[[", "rname")))
    strand <- unlist(unname(lapply(bam, "[[", "strand")))
    pos <- unlist(unname(lapply(bam, "[[", "pos")))
    cigar <-
        unlist(unname(lapply(bam, function(x) as.character(cigars(x$cigar)))))
    Alignments0(rname=rname, strand=strand, pos=pos, cigar=cigar)
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setAs("Alignments1", "Alignments0",
    function(from)
        Alignments0(rname=rname(from), strand=strand(from),
                    pos=start(from), cigar=cigar(from))

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
        if (length(i) == 0L) {
            i <- integer(0)
        } else if (is.numeric(i)) {
            if (min(i) < 0L)
                i <- seq_len(lx)[i]
            else if (!is.integer(i))
                i <- as.integer(i)
        } else if (is.logical(i)) {
            if (length(i) > lx)
                stop("subscript out of bounds")
            i <- seq_len(lx)[i]
        } else {
            stop("invalid subscript type")
        }
        x@strand <- .subsetCompactRawVector(x@strand, i)
        x@rname <- x@rname[i]
        x@cigar <- x@cigar[i]
        x@ranges <- x@ranges[i]
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "shift" method.
###

setMethod("shift", "Alignments0",
    function(x, shift, use.names=TRUE)
    {
        shift <- normargShift(shift, length(x))
        shift <- rep.int(shift, elementLengths(x@ranges))
        x@ranges@unlistData <- shift(x@ranges@unlistData,
                                     shift, use.names=use.names)
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "updateCigarAndStart" method.
###
### Performs atomic update of the cigar/start information.
###

setMethod("updateCigarAndStart", "Alignments0",
    function(x, cigar=NULL, start=NULL)
    {
        if (is.null(cigar))
            cigar <- cigar(x)
        else if (!is.character(cigar) || length(cigar) != length(x))
            stop("when not NULL, 'cigar' must be a character vector ",
                 "of the same length as 'x'")
        if (is.null(start))
            start <- start(x)
        else if (!is.integer(start) || length(start) != length(x))
            stop("when not NULL, 'start' must be an integer vector ",
                 "of the same length as 'x'")
        ranges <- cigarToIRangesListByAlignment(cigar, start)
        ## Atomic update (until the 2 slots are updated, x@cigar and x@ranges
        ## will be temporarily out of sync):
        x@cigar <- cigar
        x@ranges <- ranges
        x
    }
)

