### =========================================================================
### Alignments objects
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

setMethod("length", "Alignments0", function(x) length(x@cigar))

setMethod("rname", "Alignments0", function(x) x@rname)

setMethod("strand", "Alignments0",
    function(x)
    {
        is_minus <- .compactRawVectorAsLogical(x@strand, length(x))
        strand(ifelse(is_minus, "-", "+"))
    }
)

setMethod("cigar", "Alignments0", function(x) x@cigar)

setMethod("ranges", "Alignments0", function(x) x@ranges)

setMethod("qwidth", "Alignments0", function(x) cigarToQWidth(cigar(x)))

setMethod("isSimple", "Alignments0",
    function(x) qwidth(x) == (max(ranges(x)) - min(ranges(x)) + 1L)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.Alignments0.rname <- function(x)
{
    x_rname <- rname(x)
    if (!is.factor(x_rname) || !is.character(levels(x_rname))
     || !is.null(names(x_rname)) || any(is.na(x_rname)))
        return("'rname(x)' must be an unnamed character factor with no NAs")
    if (length(x_rname) != length(cigar(x)))
        return("'rname(x)' and 'cigar(x)' must have the same length")
    NULL
}

.valid.Alignments0.strand <- function(x)
{
    x_strand <- strand(x)
    if (!is.factor(x_strand) || !identical(levels(x_strand), levels(strand()))
     || !is.null(names(x_strand)) || any(is.na(x_strand)))
        return("'strand(x)' must be an unnamed character factor with no NAs (and with levels +, - and *)")
    if (length(x_strand) != length(cigar(x)))
        return("'strand(x)' and 'cigar(x)' must have the same length")
    NULL
}

.valid.Alignments0.cigar <- function(x)
{
    x_cigar <- cigar(x)
    if (!is.character(x_cigar) || !is.null(names(x_cigar)) || any(is.na(x_cigar)))
        return("'cigar(x)' must be an unnamed character vector with no NAs")
    tmp <- validCigar(x_cigar)
    if (!is.null(tmp))
        return(paste("in 'cigar(x)':", tmp))
    NULL
}

.valid.Alignments0.ranges <- function(x)
{
    x_ranges <- ranges(x)
    if (!is(x_ranges, "CompressedNormalIRangesList") || !is.null(names(x_ranges)))
        return("'ranges(x)' must be an unnamed CompressedNormalIRangesList object")
    if (length(x_ranges) != length(cigar(x)))
        return("'ranges(x)' and 'cigar(x)' must have the same length")
    if (any(elementLengths(x_ranges) == 0L))
        return("'ranges(x)' has elements with no ranges")
    x_ranges2 <- cigarToIRangesListByAlignment(cigar(x), min(x_ranges))
    if (!identical(x_ranges2, x_ranges))
        return("'ranges(x)' and 'cigar(x)' are incompatible")
    NULL
}

.valid.Alignments0 <- function(x)
{
    c(.valid.Alignments0.rname(x),
      .valid.Alignments0.strand(x),
      .valid.Alignments0.cigar(x),
      .valid.Alignments0.ranges(x))
}

setValidity2("Alignments0", .valid.Alignments0, where=asNamespace("Rsamtools"))


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
                          pos=min(ranges(x)),
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
                         pos=sketch(min(ranges(object))),
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
readBAMasAlignments <- function(file, index=file, which=RangesList())
{
    param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE,
                                           isDuplicate=FALSE),
                          what=c("rname", "strand", "pos", "cigar"),
                          which=which)
    bam <- scanBam(file, index=index, param=param)
    ans_rname <- unlist(unname(lapply(bam, function(x) x$rname)))
    if (!is.factor(ans_rname))
        ans_rname <- as.factor(ans_rname)
    ans_strand <- unlist(unname(lapply(bam, function(x) x$strand)))
    ans_strand <- .logicalAsCompactRawVector(ans_strand == "-")
    ans_cigar <- unlist(unname(lapply(bam, function(x) x$cigar@.cigar)))
    if (is.factor(ans_cigar))
        ans_cigar <- as.vector(ans_cigar)
    ans_pos <- unlist(unname(lapply(bam, function(x) x$pos)))
    ans_ranges <- cigarToIRangesListByAlignment(ans_cigar, ans_pos)
    new("Alignments0", rname=ans_rname, strand=ans_strand,
                       cigar=ans_cigar, ranges=ans_ranges)
}


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
### The "coverage" method.
###

setMethod("coverage", "Alignments0",
    function(x, start=NA, end=NA, shift=0L, width=NULL, weight=1L, ...)
    {
        if (!identical(start, NA) || !identical(end, NA)
         || !identical(shift, 0L) || !is.null(width) || !identical(weight, 1L))
            stop("'start', 'end', 'shift', 'width' and 'weight' arguments ",
                 "are not supported yet, sorry!")
        irl <- cigarToIRangesListByRName(cigar(x), rname(x), min(ranges(x)))
        irl <- irl[elementLengths(irl) != 0]  # drop empty elements
        coverage(irl)
    }
)
