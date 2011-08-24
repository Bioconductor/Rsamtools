scanBamFlag_test <- function()
{
    checkIdentical(formalArgs(scanBamFlag), Rsamtools:::.FLAG_BITNAMES)

    ## HP: Is this really testing something about scanBamFlag?
    isUnmappedQuery <- FALSE
    flags0 <- scanBamFlag(isUnmappedQuery=FALSE)
    flags1 <- scanBamFlag(isUnmappedQuery=isUnmappedQuery)
    checkIdentical(flags0, flags1)
}
