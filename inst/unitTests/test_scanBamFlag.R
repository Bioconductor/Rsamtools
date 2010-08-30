scanBamFlag_test <- function()
{
    isUnmappedQuery <- FALSE
    flags0 <- scanBamFlag(isUnmappedQuery=FALSE)
    flags1 <- scanBamFlag(isUnmappedQuery=isUnmappedQuery)
    checkIdentical(flags0, flags1)
}
