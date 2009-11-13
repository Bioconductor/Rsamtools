###
cigarToReadWidth <- function(cigar)
{
    if (!is.character(cigar)) {
        if (!is.factor(cigar) || !is.character(levels(cigar)))
            stop("'cigar' must be a character vector/factor")
        cigar <- as.vector(cigar)
    }
    .Call(".cigar_to_read_width", cigar, PACKAGE="Rsamtools")
}

cigarToIRanges <- function(cigar, drop.D.ranges=FALSE)
{
    if (is.factor(cigar) && is.character(levels(cigar)))
        cigar <- as.vector(cigar)
    if (!isSingleString(cigar))
        stop("'cigar' must be a single string")
    if (!isTRUEorFALSE(drop.D.ranges))
        stop("'drop.D.ranges' must be TRUE or FALSE")
    .Call(".cigar_to_IRanges", cigar, drop.D.ranges, PACKAGE="Rsamtools")
}

### NOTE: 'strand' is ignored for now.
cigarToIRangesList <- function(cigar, rname, strand, pos, drop.D.ranges=FALSE)
{
    if (!is.character(cigar)) {
        if (!is.factor(cigar) || !is.character(levels(cigar)))
            stop("'cigar' must be a character vector/factor")
        cigar <- as.vector(cigar)
    }
    if (!is.factor(rname) || !is.character(levels(rname))) {
        if (!is.character(rname))
            stop("'rname' must be a character vector/factor")
        rname <- as.factor(rname)
    }
    if (!is.numeric(pos))
        stop("'pos' must be a vector of integers")
    if (!is.integer(pos))
        pos <- as.integer(pos)
    if (length(cigar) != length(rname) || length(rname) != length(pos))
        stop("'cigar', 'rname' and 'pos' must have the same length")
    if (!isTRUEorFALSE(drop.D.ranges))
        stop("'drop.D.ranges' must be TRUE or FALSE")
    C_ans <- .Call(".cigar_to_list_of_IRanges",
                   cigar, rname, strand, pos, drop.D.ranges,
                   PACKAGE="Rsamtools")
    if (length(C_ans) < 100L)
        IRangesList(C_ans, compress=FALSE)
    else
        IRangesList(C_ans, compress=TRUE)
}

