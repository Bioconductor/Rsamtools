### =========================================================================
### GappedAlignments objects
### -------------------------------------------------------------------------
###
### Formal API:
###   length(x)   - single integer. Nb of alignments in 'x'.
###   rname(x)    - character factor of the same length as 'x'.
###   rname(x) <- value - replacement form of 'rname(x)'.
###   strand(x)   - character factor of the same length as 'x' (levels: +, -,
###                 *).
###   cigar(x)    - character vector of the same length as 'x'.
###   qwidth(x)   - integer vector of the same length as 'x'.
###   granges(x)  - GRangesList object of the same length as 'x'.
###   ranges(x)   - CompressedNormalIRangesList object of the same length as
###                 'x'.
###   start(x), end(x), width(x) - integer vectors of the same length as 'x'.
###   ngap(x)     - integer vector of the same length as 'x'.
###   as.data.frame(x) - just a convenience used by show(x).
###   show(x)     - compact display in a data.frame-like fashion.
###   x[i]        - GappedAlignments object of the same class as 'x'
###                 (endomorphism).
###
###   shift(x, shift) - GappedAlignments object of the same length and class
###                 as 'x' (endomorphism).
###
###   updateCigarAndStart(x, cigar=NULL, start=NULL) - GappedAlignments
###                 object of the same length and class as 'x' (endomorphism).
###
###   qnarrow(x, start=NA, end=NA, width=NA) - GappedAlignments object of the
###                 same length and class as 'x' (endomorphism).
###
###   narrow(x, start=NA, end=NA, width=NA) - GappedAlignments object of the
###                 same length and class as 'x' (endomorphism).
###
###   coverage(x) - named RleList object with one element (integer-Rle) per
###                 unique reference sequence.
###
###   findOverlaps(query, subject) - 'query' or 'subject' or both are
###                 GappedAlignments objects. Just a convenient wrapper for
###                 'findOverlaps(granges(query), subject, ...)', etc...
###
### Concrete GappedAlignments implementations just need to implement:
###   length, rname, rname<-, strand, cigar, granges, ranges,
###   [, shift and updateCigarAndStart
### and the default methods defined in this file will work.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor-like methods.
###

setMethod("length", "GappedAlignments", function(x) length(x@cigar))

setMethod("cigar", "GappedAlignments", function(x) x@cigar)
setMethod("qwidth", "GappedAlignments", function(x) cigarToQWidth(x@cigar))

setMethod("start", "GappedAlignments", function(x, ...) min(ranges(x)))
setMethod("end", "GappedAlignments", function(x, ...) max(ranges(x)))
setMethod("width", "GappedAlignments", function(x) cigarToWidth(x@cigar))

setMethod("ngap", "GappedAlignments",
    function(x) {elementLengths(ranges(x)) - 1L}
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Used by the "rname<-" methods.
###

normargRNameReplaceValue <- function(x, value, ans.type=c("factor", "Rle"))
{
    ans.type <- match.arg(ans.type)
    if (!is.factor(value)
     && !is.character(value)
     && (!is(value, "Rle") || !is.character(runValue(value))
                              && !is.factor(runValue(value))))
        stop("'rname' value must be a character factor/vector, ",
             "or a 'character' Rle, or a 'factor' Rle")
    if (ans.type == "factor" && !is.factor(value)) {
        value <- as.factor(value)
    } else if (ans.type == "Rle") {
        ## We want to return a 'character' Rle.
        if (!is(value, "Rle")) {
            if (!is.character(value))
                value <- as.character(value)
            value <- Rle(value)
        } else if (!is.character(runValue(value))) {
            runValue(value) <- as.character(runValue(value))
        }
    }
    if (length(value) != length(x))
        stop("'rname' value must be the same length as the object")

    ## Check the mapping between old and new values.
    old <- rname(x)  # guaranteed to be factor
    if (ans.type == "factor") {
        new <- value
    } else {
        old <- Rle(old)
        if (!identical(runLength(old), runLength(value))) {
            warning("mapping between old an new 'rname' values ",
                    "is not one-to-one")
            return(value)
        }
        old <- runValue(old)
        new <- runValue(value)
    }
    mapping <- unique(data.frame(old=old, new=new))
    if (any(duplicated(mapping$old)) || any(duplicated(mapping$new)))
        warning("mapping between old an new 'rname' values ",
                "is not one-to-one")
    value
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper function.
###

### The arguments are the different components of a GappedAlignments object
### instead of just the GappedAlignments object itself (arg 'x'). This allows
### the function to be used in different contexts e.g. within the Alignments1()
### constructor when 'x' doesn't exist yet but is in the process of being
### constructed.
GappedAlignmentsAsGRangesList <- function(rname, strand, ranges, check=TRUE)
{
    if (check) {
        if (is.factor(rname) && is.character(levels(rname)))
            rname <- as.character(rname)
        else if (!is.character(rname))
            stop("'rname' must be a character vector/factor")
        if (any(is.na(rname)))
            stop("'rname' cannot have NAs")
        if (!is.factor(strand) || !identical(levels(strand), .STRAND_LEVELS))
            stop("'strand' must be a character factor")
        if (!is(ranges, "CompressedNormalIRangesList"))
            stop("'ranges' must be a CompressedNormalIRangesList object")
    }
    nrg_per_alignment <- elementLengths(ranges)
    seqnames <- Rle(rname, nrg_per_alignment)
    strand <- Rle(strand, nrg_per_alignment)
    granges <- GRanges(seqnames=seqnames,
                       ranges=ranges@unlistData,
                       strand=strand)
    new("GRangesList", unlistData=granges, partitioning=ranges@partitioning)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.GappedAlignments.cigar <- function(x)
{
    x_cigar <- cigar(x)
    if (!is.character(x_cigar) || !is.null(names(x_cigar)) || any(is.na(x_cigar)))
        return("'cigar(x)' must be an unnamed character vector with no NAs")
    tmp <- validCigar(x_cigar)
    if (!is.null(tmp))
        return(paste("in 'cigar(x)':", tmp))
    NULL
}

.valid.GappedAlignments.rname <- function(x)
{
    x_rname <- rname(x)
    if (!is.factor(x_rname) || !is.character(levels(x_rname))
     || !is.null(names(x_rname)) || any(is.na(x_rname)))
        return("'rname(x)' must be an unnamed character factor with no NAs")
    if (length(x_rname) != length(cigar(x)))
        return("'rname(x)' and 'cigar(x)' must have the same length")
    NULL
}

.valid.GappedAlignments.strand <- function(x)
{
    x_strand <- strand(x)
    if (!is.factor(x_strand) || !identical(levels(x_strand), levels(strand()))
     || !is.null(names(x_strand)) || any(is.na(x_strand)))
        return("'strand(x)' must be an unnamed character factor with no NAs (and with levels +, - and *)")
    if (length(x_strand) != length(cigar(x)))
        return("'strand(x)' and 'cigar(x)' must have the same length")
    NULL
}

.valid.GappedAlignments.ranges <- function(x)
{
    x_ranges <- ranges(x)
    if (!is(x_ranges, "CompressedNormalIRangesList") || !is.null(names(x_ranges)))
        return("'ranges(x)' must be an unnamed CompressedNormalIRangesList object")
    if (length(x_ranges) != length(cigar(x)))
        return("'ranges(x)' and 'cigar(x)' must have the same length")
    if (any(elementLengths(x_ranges) == 0L))
        return("'ranges(x)' has elements with no ranges")
    x_start <- min(x_ranges)
    x_end <- max(x_ranges)
    if (!identical(x_end - x_start + 1L, width(x)))
        return("'ranges(x)' and 'width(x)' are out of sync")
    x_ranges2 <- cigarToIRangesListByAlignment(cigar(x), x_start)
    if (!identical(x_ranges2, x_ranges))
        return("'ranges(x)' and 'cigar(x)' are out of sync")
    NULL
}

.valid.GappedAlignments <- function(x)
{
    c(.valid.GappedAlignments.cigar(x),
      .valid.GappedAlignments.rname(x),
      .valid.GappedAlignments.strand(x),
      #.valid.GappedAlignments.granges(x),
      .valid.GappedAlignments.ranges(x))
}

setValidity2("GappedAlignments", .valid.GappedAlignments,
             where=asNamespace("Rsamtools"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "as.data.frame" and "show" methods.
###

setMethod("as.data.frame", "GappedAlignments",
    function(x, row.names=NULL, optional=FALSE, ...)
    {
        if (!(is.null(row.names) || is.character(row.names)))
            stop("'row.names' must be NULL or a character vector")
        ans <- data.frame(rname=rname(x),
                          strand=strand(x),
                          cigar=cigar(x),
                          qwidth=qwidth(x),
                          start=start(x),
                          end=end(x),
                          width=width(x),
                          ngap=ngap(x),
                          row.names=row.names,
                          check.rows=TRUE,
                          check.names=FALSE,
                          stringsAsFactors=FALSE)
        return(ans)
    }
)

setMethod("show", "GappedAlignments",
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
                         qwidth=sketch(qwidth(object)),
                         start=sketch(start(object)),
                         end=sketch(end(object)),
                         width=sketch(width(object)),
                         ngap=sketch(ngap(object)),
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
### Subsetting.
###

### Doesn't subset anything, just convert 'i' into integer vector and return
### it. Supported input types for 'i': numeric vector, logical vector, NULL
### and missing.
setMethod("[", "GappedAlignments",
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
        i
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "qnarrow" and "narrow" methods.
###

setMethod("qnarrow", "GappedAlignments",
    function(x, start=NA, end=NA, width=NA)
    {
        ans_cigar <- cigarQNarrow(cigar(x),
                                  start=start, end=end, width=width)
        ans_start <- start(x) + attr(ans_cigar, "rshift")
        updateCigarAndStart(x, cigar=ans_cigar, start=ans_start)
    }
)

setMethod("narrow", "GappedAlignments",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        ans_cigar <- cigarNarrow(cigar(x),
                                 start=start, end=end, width=width)
        ans_start <- start(x) + attr(ans_cigar, "rshift")
        updateCigarAndStart(x, cigar=ans_cigar, start=ans_start)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "coverage" method.
###

setMethod("coverage", "GappedAlignments",
    function(x, start=NA, end=NA, shift=0L, width=NULL, weight=1L, ...)
    {
        if (!identical(start, NA) || !identical(end, NA)
         || !identical(shift, 0L) || !is.null(width) || !identical(weight, 1L))
            stop("'start', 'end', 'shift', 'width' and 'weight' arguments ",
                 "are not supported yet, sorry!")
        irl <- cigarToIRangesListByRName(cigar(x), rname(x), start(x))
        irl <- irl[elementLengths(irl) != 0]  # drop empty elements
        coverage(irl)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "findOverlaps" methods.
###

setMethod("findOverlaps", c("GappedAlignments", "ANY"),
    function(query, subject, maxgap=0, multiple=TRUE,
             type=c("any", "start", "end"))
    {
        callGeneric(granges(query), subject,
                    maxgap=maxgap, multiple=multiple, type=type)
    }
)

setMethod("findOverlaps", c("ANY", "GappedAlignments"),
    function(query, subject, maxgap=0, multiple=TRUE,
             type=c("any", "start", "end"))
    {
        callGeneric(query, granges(subject),
                    maxgap=maxgap, multiple=multiple, type=type)
    }
)

### Not strictly needed! Defining the above 2 methods covers that case but
### with the following note:
###   > findOverlaps(al1, al0)
###   Note: Method with signature "GappedAlignments#ANY" chosen for
###    function "findOverlaps", target signature "Alignments1#Alignments0".
###    "ANY#GappedAlignments" would also be valid
setMethod("findOverlaps", c("GappedAlignments", "GappedAlignments"),
    function(query, subject, maxgap=0, multiple=TRUE,
             type=c("any", "start", "end"))
    {
        callGeneric(granges(query), granges(subject),
                    maxgap=maxgap, multiple=multiple, type=type)
    }
)

