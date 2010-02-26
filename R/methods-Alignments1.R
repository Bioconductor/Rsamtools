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
        xgrg <- x@grglist
        as.factor(seqnames(xgrg@unlistData))[xgrg@partitioning@end]
    }
)

setReplaceMethod("rname", "Alignments1",
    function(x, value)
    {
        value <- normargRNameReplaceValue(x, value, ans.type="Rle")
        value <- rep.int(value, elementLengths(x@grglist))
        seqnames(x@grglist@unlistData) <- value
        x
    }
)

setMethod("strand", "Alignments1",
    function(x)
    {
        xgrg <- x@grglist
        as.factor(strand(xgrg@unlistData))[xgrg@partitioning@end]
    }
)

setMethod("grglist", "Alignments1", function(x) x@grglist)

setMethod("rglist", "Alignments1",
    function(x) as(ranges(x@grglist), "CompressedNormalIRangesList")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors.
###

Alignments1 <- function(rname=factor(), strand=BSgenome::strand(),
                        pos=integer(), cigar=character())
{
    rglist <- cigarToIRangesListByAlignment(cigar, pos)
    grglist <- GappedAlignmentsAsGRangesList(rname, strand, rglist)
    new("Alignments1", cigar=cigar, grglist=grglist)
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

setAs("GappedAlignments", "Alignments1",
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
        i <- callNextMethod()
        x@cigar <- x@cigar[i]
        x@grglist <- x@grglist[i]
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
        rglist <- cigarToIRangesListByAlignment(cigar, start)
        grglist <- GappedAlignmentsAsGRangesList(rname(x), strand(x), rglist)
        ## Atomic update (until the 2 slots are updated, x@cigar and x@grglist
        ## will be temporarily out of sync):
        x@cigar <- cigar
        x@grglist <- grglist
        x
    }
)

