### =========================================================================
### Alignments1 objects
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor-like methods.
###

setMethod("rname", "Alignments1",
    function(x)
    {
        xgrg <- x@granges
        as.factor(seqnames(xgrg@unlistData))[xgrg@partitioning@end]
    }
)

setMethod("strand", "Alignments1",
    function(x)
    {
        xgrg <- x@granges
        as.factor(strand(xgrg@unlistData))[xgrg@partitioning@end]
    }
)

setMethod("granges", "Alignments1", function(x) x@granges)

setMethod("ranges", "Alignments1",
    function(x) as(ranges(x@granges), "CompressedNormalIRangesList")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.Alignments1.cigar <- function(x)
{
    x_cigar <- cigar(x)
    if (!is.character(x_cigar) || !is.null(names(x_cigar)) || any(is.na(x_cigar)))
        return("'cigar(x)' must be an unnamed character vector with no NAs")
    tmp <- validCigar(x_cigar)
    if (!is.null(tmp))
        return(paste("in 'cigar(x)':", tmp))
    NULL
}

.valid.Alignments1.ranges <- function(x)
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

.valid.Alignments1 <- function(x)
{
    c(#.valid.Alignments1.rname(x),
      #.valid.Alignments1.strand(x),
      .valid.Alignments1.cigar(x))
}

setValidity2("Alignments1", .valid.Alignments1,
    where=asNamespace("Rsamtools"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors.
###

Alignments1 <- function(rname=factor(), strand=BSgenome::strand(),
                        pos=integer(), cigar=character())
{
    ranges <- cigarToIRangesListByAlignment(cigar, pos)
    granges <- GappedAlignmentsAsGRangesList(rname, strand, ranges)
    new("Alignments1", cigar=cigar, granges=granges)
}

setMethod(readBAMasAlignments1, "character",
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
        return(Alignments1())
    rname <- unlist(unname(lapply(bam, "[[", "rname")))
    strand <- unlist(unname(lapply(bam, "[[", "strand")))
    pos <- unlist(unname(lapply(bam, "[[", "pos")))
    cigar <-
        unlist(unname(lapply(bam, function(x) as.character(cigars(x$cigar)))))
    Alignments1(rname=rname, strand=strand, pos=pos, cigar=cigar)
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setAs("Alignments0", "Alignments1",
    function(from)
        Alignments1(rname=rname(from), strand=strand(from),
                    pos=start(from), cigar=cigar(from))

)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

### Supported 'i' types: numeric vector, logical vector, NULL and missing.
setMethod("[", "Alignments1",
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
        x@cigar <- x@cigar[i]
        x@granges <- x@granges[i]
        x
    }
)

