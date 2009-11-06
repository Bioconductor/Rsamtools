### NOTE: 'strand' is ignored for now.
cigarToIRangesList <- function(rname, strand, pos, cigar)
{
    if (is.factor(rname))
        rname <- as.vector(rname)
    if (!is.character(rname))
        stop("'rname' must be a character vector (or a factor)")
    if (!is.numeric(pos))
        stop("'pos' must be a vector of integers")
    if (!is.integer(pos))
        pos <- as.integer(pos)
    if (is.factor(cigar))
        cigar <- as.vector(cigar)
    if (!is.character(cigar))
        stop("'cigar' must be a character vector (or a factor)")
    rname2rank_env <- new.env(hash=TRUE, parent=emptyenv())
    C_ans <- .Call(".cigar_to_list_of_IRanges",
                   rname, strand, pos, cigar, rname2rank_env,
                   PACKAGE="Rsamtools")
    IRangesList(C_ans, compress=TRUE)
}

