### =========================================================================
### Alignments2 objects
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor-like methods.
###

setMethod("rname", "Alignments2", function(x) x@rname)

setReplaceMethod("rname", "Alignments2",
    function(x, value)
    {
        x@rname <- normargRNameReplaceValue(x, value)
        x
    }
)

setMethod("strand", "Alignments2",
    function(x) strand(.compactRawVectorAsLogical(x@strand, length(x)))
)

setMethod("granges", "Alignments2",
    function(x)
        GappedAlignmentsAsGRangesList(x@rname, strand(x), ranges(x))
)

setMethod("ranges", "Alignments2",
    function(x)
        cigarToIRangesListByAlignment(x@cigar, x@start)
)

setMethod("start", "Alignments2", function(x, ...) x@start)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors.
###

Alignments2 <- function(rname=factor(), strand=BSgenome::strand(),
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
    if (!is.integer(pos) || any(is.na(pos)))
        stop("'pos' must be an integer vector with no NAs")
    if (!is.character(cigar) || any(is.na(cigar)))
        stop("'cigar' must be a character vector with no NAs")
    strand <- .logicalAsCompactRawVector(strand == "-")
    new("Alignments2", rname=rname, strand=strand,
                       cigar=cigar, start=pos)
}

### This is our only constructor for now.
setMethod(readBAMasAlignments2, "character", 
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
        return(Alignments2())
    rname <- unlist(unname(lapply(bam, "[[", "rname")))
    strand <- unlist(unname(lapply(bam, "[[", "strand")))
    pos <- unlist(unname(lapply(bam, "[[", "pos")))
    cigar <-
        unlist(unname(lapply(bam, function(x) as.character(cigars(x$cigar)))))
    Alignments2(rname=rname, strand=strand, pos=pos, cigar=cigar)
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setAs("GappedAlignments", "Alignments2",
    function(from)
        Alignments2(rname=rname(from), strand=strand(from),
                    pos=start(from), cigar=cigar(from))

)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

### Supported 'i' types: numeric vector, logical vector, NULL and missing.
setMethod("[", "Alignments2",
    function(x, i, j, ... , drop=TRUE)
    {
        i <- callNextMethod()
        x@cigar <- x@cigar[i]
        x@rname <- x@rname[i]
        x@strand <- .subsetCompactRawVector(x@strand, i)
        x@start <- x@start[i]
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "updateCigarAndStart" method.
###
### Performs atomic update of the cigar/start information.
###

setMethod("updateCigarAndStart", "Alignments2",
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
        x@cigar <- cigar
        x@start <- start
        x
    }
)

