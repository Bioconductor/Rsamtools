### =========================================================================
### GappedAlignments objects
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor-like methods.
###

setMethod("rname", "GappedAlignments",
    function(x) as.factor(seqnames(x@unlistData)[x@partitioning])
)

setMethod("strand", "GappedAlignments",
    function(x) as.factor(strand(x@unlistData)[x@partitioning])
)

setMethod("cigar", "GappedAlignments", function(x) x@cigar)

setMethod("ranges", "GappedAlignments",
    function(x) as(callNextMethod(), "CompressedNormalIRangesList")
)

setMethod("qwidth", "GappedAlignments", function(x) cigarToQWidth(cigar(x)))


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

.valid.GappedAlignments.ranges <- function(x)
{
    x_ranges <- ranges(x)
    if (length(x_ranges) != length(cigar(x)))
        return("'ranges(x)' and 'cigar(x)' must have the same length")
    if (any(elementLengths(x_ranges) == 0L))
        return("'ranges(x)' has elements with no ranges")
    x_ranges2 <- cigarToIRangesListByAlignment(cigar(x), min(x_ranges))
    if (!identical(x_ranges2, x_ranges))
        return("'ranges(x)' and 'cigar(x)' are incompatible")
    NULL
}

.valid.GappedAlignments <- function(x)
{
    c(#.valid.GappedAlignments.rname(x),
      #.valid.GappedAlignments.strand(x),
      .valid.GappedAlignments.cigar(x))
}

setValidity2("GappedAlignments", .valid.GappedAlignments,
    where=asNamespace("Rsamtools"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

setMethod("as.data.frame", "GappedAlignments",
    function(x, row.names=NULL, optional=FALSE, ...)
    {
        if (!(is.null(row.names) || is.character(row.names)))
            stop("'row.names' must be NULL or a character vector")
        ans <- data.frame(#rname=rname(x),
                          #strand=strand(x),
                          cigar=cigar(x),
                          pos=min(ranges(x)),
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

### This helper constructor is not exported for now.
.GappedAlignments <- function(rname=factor(), strand=strand(),
                              pos=integer(0), cigar=character(0))
{
    if (is.factor(rname) && is.character(levels(rname)))
        rname <- as.character(rname)
    else if (!is.character(rname)) 
        stop("'rname' must be a character vector/factor")
    if (any(is.na(rname)))
        stop("'rname' cannot have NAs")
    if (!is.factor(strand) || !is.character(levels(strand)))
        stop("'strand' must be a character factor")
    if (!is.character(cigar) || any(is.na(cigar)))
        stop("'cigar' must be a character vector with no NAs")
    rg_list <- cigarToIRangesListByAlignment(cigar, pos)
    nrg_per_alignment <- elementLengths(rg_list)
    rname <- Rle(rname, nrg_per_alignment)
    strand <- Rle(strand, nrg_per_alignment)
    ranges <- unlist(rg_list)
    unlistData <- GenomicFeature(seqnames=rname, ranges=ranges, strand=strand)
    partitioning <- PartitioningByEnd(cumsum(nrg_per_alignment))
    new("GappedAlignments",
        unlistData=unlistData,
        partitioning=partitioning,
        cigar=cigar)
}

setMethod(readBAMasGappedAlignments, "character",
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
        return(.GappedAlignments())
    rname <- unlist(unname(lapply(bam, "[[", "rname")))
    strand <- unlist(unname(lapply(bam, "[[", "strand")))
    pos <- unlist(unname(lapply(bam, "[[", "pos")))
    cigar <-
        unlist(unname(lapply(bam, function(x) as.character(cigars(x$cigar)))))
    .GappedAlignments(rname=rname, strand=strand, pos=pos, cigar=cigar)
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

### Supported 'i' types: numeric vector, logical vector, NULL and missing.
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
        x <- callNextMethod()
        x@cigar <- x@cigar[i]
        x
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
        irl <- cigarToIRangesListByRName(cigar(x), rname(x), min(ranges(x)))
        irl <- irl[elementLengths(irl) != 0]  # drop empty elements
        coverage(irl)
    }
)
