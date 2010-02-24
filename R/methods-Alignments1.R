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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "shift" method.
###

setMethod("shift", "Alignments1",
    function(x, shift, use.names=TRUE)
    {
        shift <- normargShift(shift, length(x))
        shift <- rep.int(shift, elementLengths(x@granges))
        x@granges@unlistData <- shift(x@granges@unlistData,
                                      shift, use.names=use.names)
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "updateCigarAndStart" method.
###
### Performs atomic update of the cigar/start information.
###

setMethod("updateCigarAndStart", "Alignments1",
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
        granges <- GappedAlignmentsAsGRangesList(rname(x), strand(x), ranges)
        ## Atomic update (until the 2 slots are updated, x@cigar and x@granges
        ## will be temporarily out of sync):
        x@cigar <- cigar
        x@granges <- granges
        x
    }
)

