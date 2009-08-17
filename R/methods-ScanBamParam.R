ScanBamParam <-
    function(flag=scanBamFlag(), simpleCigar=FALSE,
             what=scanBamWhat(), which=RangesList())
{
    new("ScanBamParam", flag=flag, simpleCigar=simpleCigar,
        what=what, which=which)
}

setValidity("ScanBamParam", function(object) {
    msg <- NULL
    flag <- bamFlag(object)
    simpleCigar <- bamSimpleCigar(object)
    which <- bamWhich(object)
    what <- bamWhat(object)
    if (length(flag) != 2 || typeof(flag) != "integer")
        msg <- c(msg, "'flag' must be integer(2)")
    else {
        if (!identical(names(flag), c("keep0", "keep1")))
            msg <- c(msg, "'names(flag)' must be c('keep0', 'keep1')")
        else if (any(flag < 0 | flag >= 2^16))
            msg <- c(msg, "'flag' values must be >=0, <2048")
    }
    if (length(simpleCigar) != 1 && !is.na(simpleCigar))
        msg <- c(msg, "'simpleCigar' must be logical(1), not NA")
    if (length(which) != 0 &&
        (any(!nzchar(names(which))) || any(is.na(names(which)))))
        msg <- c(msg, "'which' elements must be named (not NA)")
    if (!all(what %in% scanBamWhat()))
        msg <- c(msg, "all 'what' must be from 'scanBamWhat()'")
    if (is.null(msg)) TRUE else msg
})

bamFlag <- function(object) slot(object, "flag")
bamSimpleCigar <- function(object) slot(object, "simpleCigar")
bamWhich <- function(object) slot(object, "which")
bamWhat <- function(object) slot(object, "what")

## helpers 

scanBamFlag <-
    function(isPaired=NA, isProperPair=NA, isUnmappedQuery=NA,
             hasUnmappedMate=NA, isMinusStrand=NA, isMateStrand=NA,
             isFirstMateRead=NA, isSecondMateRead=NA, # redundant
             isPrimaryRead=NA, isValidVendorRead=NA, isDuplicate=NA)
    ## NA: keep either 0 or 1 flag; FALSE: keep 0 flag; TRUE: keep 1 flag
{
    flag <- c(isPaired=1L, isProperPair=2L, isUnmappedQuery=4L,
              hasUnmappedMate=8L, isMinusStrand=16L, isMateStrand=32L,
              isFirstMateRead=64L, isSecondMateRead=128L,
              isPrimaryRead=256L, isValidVendorRead=512L, isDuplicate=1024L)
    args <- as.list(match.call())[-1]
    if (any(sapply(args, length) > 1L))
        stop("all arguments must be logical(1)")
    if (length(args) > 0) {
        ## keep0: NA | FALSE --> drop !NA & TRUE
        idx <- names(args[sapply(args, function(x) !is.na(x) && x)])
        keep0 <- Reduce("+", flag[ !names(flag) %in% idx ], 0L)
        ## keep1: NA | TRUE --> drop !NA & FALSE
        idx <- names(args[sapply(args, function(x) !is.na(x) && !x)])
        keep1 <- Reduce("+", flag[ !names(flag) %in% idx ], 0L)
    } else
        keep0 <- keep1 <- Reduce("+", flag, 0L)
    c(keep0=keep0, keep1=keep1)
}

scanBamWhat <- function()
{
    names(.Call(.scan_bam_template))
}

