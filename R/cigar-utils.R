###
cigarToReadWidth <- function(cigar, after.hard.clipping=FALSE)
{
    if (!is.character(cigar)) {
        if (!is.factor(cigar) || !is.character(levels(cigar)))
            stop("'cigar' must be a character vector/factor")
        cigar <- as.vector(cigar)
    }
    if (!isTRUEorFALSE(after.hard.clipping))
        stop("'after.hard.clipping' must be TRUE or FALSE")
    .Call(".cigar_to_read_width",
          cigar, after.hard.clipping,
          PACKAGE="Rsamtools")
}

cigarToIRanges <- function(cigar, drop.D.ranges=FALSE, merge.ranges=TRUE)
{
    if (is.factor(cigar) && is.character(levels(cigar)))
        cigar <- as.vector(cigar)
    if (!isSingleString(cigar))
        stop("'cigar' must be a single string")
    if (!isTRUEorFALSE(drop.D.ranges))
        stop("'drop.D.ranges' must be TRUE or FALSE")
    if (!isTRUEorFALSE(merge.ranges))
        stop("'merge.ranges' must be TRUE or FALSE")
    .Call(".cigar_to_IRanges", cigar, drop.D.ranges, merge.ranges,
          PACKAGE="Rsamtools")
}

cigarToIRangesList <- function(cigar, rname, pos, flag=NA,
                               drop.D.ranges=FALSE, merge.ranges=TRUE)
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
    if (length(cigar) != length(rname) || length(cigar) != length(pos))
        stop("'cigar', 'rname' and 'pos' must have the same length")
    if (is.vector(flag) && length(flag) == 1 && is.na(flag)) {
        ## 'flag' is a single NA of any type
        flag <- NULL
    } else {
        if (!is.numeric(flag))
            stop("'flag' must be NA or a vector of integers")
        if (!is.integer(flag))
            flag <- as.integer(flag)
        if (length(cigar) != length(flag))
            stop("'cigar' and 'flag' must have the same length")
    }
    if (!isTRUEorFALSE(drop.D.ranges))
        stop("'drop.D.ranges' must be TRUE or FALSE")
    if (!isTRUEorFALSE(merge.ranges))
        stop("'merge.ranges' must be TRUE or FALSE")
    C_ans <- .Call(".cigar_to_list_of_IRanges",
                   cigar, rname, pos, flag, drop.D.ranges, merge.ranges,
                   PACKAGE="Rsamtools")
    if (length(C_ans) < 200L)
        IRangesList(C_ans, compress=FALSE)
    else
        IRangesList(C_ans, compress=TRUE)
}

