### 
cigarToIRanges <- function(cigar)
{
    if (is.factor(cigar) && is.character(levels(cigar)))
        cigar <- as.vector(cigar)
    if (!isSingleString(cigar))
        stop("'cigar' must be a single string")
    .Call(".cigar_to_IRanges", cigar, PACKAGE="Rsamtools")
}

### NOTE: 'strand' is ignored for now.
cigarToIRangesList <- function(cigar, rname, strand, pos)
{
    if (is.factor(cigar))
        cigar <- as.vector(cigar)
    if (!is.character(cigar))
        stop("'cigar' must be a character vector/factor")
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
    C_ans <- .Call(".cigar_to_list_of_IRanges",
                   cigar, rname, strand, pos,
                   PACKAGE="Rsamtools")
    if (length(C_ans) < 100L)
        IRangesList(C_ans, compress=FALSE)
    else
        IRangesList(C_ans, compress=TRUE)
}

