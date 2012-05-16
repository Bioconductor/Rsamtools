setMethod(ScanBamParam, c(which="missing"),
          function(flag=scanBamFlag(), simpleCigar=FALSE,
                   reverseComplement=FALSE, tag=character(0),
                   what=character(0), which)
{
    which <- IRangesList()
    names(which) <- character()
    new("ScanBamParam", flag=flag, simpleCigar=simpleCigar,
        reverseComplement=reverseComplement, tag=tag, what=what,
        which=which)
})

setMethod(ScanBamParam, c(which="RangesList"),
          function(flag=scanBamFlag(), simpleCigar=FALSE,
                   reverseComplement=FALSE, tag=character(0),
                   what=character(0), which)
{
    new("ScanBamParam", flag=flag, simpleCigar=simpleCigar,
        reverseComplement=reverseComplement, tag=tag, what=what,
        which=which)
})

setMethod(ScanBamParam, c(which="RangedData"),
          function(flag=scanBamFlag(), simpleCigar=FALSE,
                   reverseComplement=FALSE, tag=character(0),
                   what=character(0), which)
{
    callGeneric(flag=flag, simpleCigar=simpleCigar,
                reverseComplement=reverseComplement, tag=tag, what=what,
                which=ranges(which))
})

setMethod(ScanBamParam, c(which="GRanges"),
          function(flag=scanBamFlag(), simpleCigar=FALSE,
                   reverseComplement=FALSE, tag=character(0),
                   what=character(0), which)
{
    callGeneric(flag=flag, simpleCigar=simpleCigar,
                reverseComplement=reverseComplement, tag=tag, what=what,
                which=split(ranges(which), seqnames(which)))
})

setValidity("ScanBamParam", function(object) {
    msg <- NULL
    flag <- bamFlag(object, asInteger=TRUE)
    simpleCigar <- bamSimpleCigar(object)
    reverseComplement <- bamReverseComplement(object)
    tag <- bamTag(object)
    what <- bamWhat(object)
    which <- bamWhich(object)
    if (length(flag) != 2 || typeof(flag) != "integer")
        msg <- c(msg, "'flag' must be integer(2)")
    else {
        if (!identical(names(flag), c("keep0", "keep1")))
            msg <- c(msg, "'names(flag)' must be c('keep0', 'keep1')")
        else if (any(flag < 0 | flag >= 2^16))
            msg <- c(msg, "'flag' values must be >=0, <2048")
    }
    if (!((1L == length(simpleCigar)) &&
          !is.na(simpleCigar)))
        msg <- c(msg, "'simpleCigar' must be logical(1), not NA")
    if (!((1L == length(reverseComplement)) &&
          !is.na(reverseComplement)))
        msg <- c(msg, "'reverseComplement' must be logical(1)")
    if (0 != length(tag) && any(2 != nchar(tag)))
        msg <- c(msg, "'tag' must be two letters, e.g,, 'MD'")
    if (!all(what %in% scanBamWhat()))
        msg <- c(msg, "'what' must be from 'scanBamWhat()'")
    if (length(which) != 0 &&
        (any(!nzchar(names(which))) || any(is.na(names(which)))))
        msg <- c(msg, "'which' elements must be named (not NA)")
    if (is.null(msg)) TRUE else msg
})

bamFlag <- function(object, asInteger=FALSE) {
    if (asInteger) {
        slot(object, "flag")
    } else {
        keep <- slot(object, "flag")
        keep0 <- bamFlagAsBitMatrix(keep[[1]])[1L,] == 1L
        flag <- keep1 <- bamFlagAsBitMatrix(keep[[2]])[1L,] == 1L
        flag[keep1 & sapply(keep0, isTRUE)] <- NA
        flag
    }
}

"bamFlag<-" <- function(object, value)
    ## FIXME: make this support input like bamFlag
{
    slot(object, "flag") <- value
    object
}
bamSimpleCigar <- function(object) slot(object, "simpleCigar")
"bamSimpleCigar<-" <- function(object, value)
{
    slot(object, "simpleCigar") <- value
    object
}
bamReverseComplement <-
    function(object) slot(object, "reverseComplement")
"bamReverseComplement<-" <- function(object, value)
{
    slot(object, "reverseComplement") <- value
    object
}
bamTag <- function(object) slot(object, "tag")
"bamTag<-" <- function(object, value)
{
    slot(object, "tag") <- value
    object
}

bamWhich <- function(object) slot(object, "which")

setGeneric("bamWhich<-",
           function(object, value) standardGeneric("bamWhich<-"))
        
setReplaceMethod("bamWhich", c("ScanBamParam", "RangesList"),
    function(object, value) 
{
    slot(object, "which") <- value
    object
})

setReplaceMethod("bamWhich", c("ScanBamParam", "GRanges"),
    function(object, value) 
{
    callGeneric(object, split(ranges(value), seqnames(value)))
})
        
setReplaceMethod("bamWhich", c("ScanBamParam", "RangedData"),
    function(object, value) 
{
    callGeneric(object, ranges(value))
})
        
setReplaceMethod("bamWhich", c("ScanBamParam", "ANY"),
    function(object, value) 
{
    callGeneric(object, as(value, "RangesList"))
})

bamWhat <- function(object) slot(object, "what")
"bamWhat<-" <- function(object, value)
{
    slot(object, "what") <- value
    object
}

## helpers 

FLAG_BITNAMES <- c(
    "isPaired",
    "isProperPair",
    "isUnmappedQuery",
    "hasUnmappedMate",
    "isMinusStrand",
    "isMateMinusStrand",
    "isFirstMateRead",
    "isSecondMateRead",
    "isNotPrimaryRead",
    "isNotPassingQualityControls",
    "isDuplicate"
)

scanBamFlag <-
    function(isPaired=NA, isProperPair=NA, isUnmappedQuery=NA,
             hasUnmappedMate=NA, isMinusStrand=NA,
             isMateMinusStrand=NA, isFirstMateRead=NA,
             isSecondMateRead=NA, # redundant
             isNotPrimaryRead=NA, isNotPassingQualityControls=NA,
             isDuplicate=NA, isValidVendorRead=NA)

    ## NA: keep either 0 or 1 flag; FALSE: keep 0 flag; TRUE: keep 1 flag
{
    flag <- IRanges:::makePowersOfTwo(length(FLAG_BITNAMES))
    names(flag) <- FLAG_BITNAMES
    args <- lapply(as.list(match.call())[-1], eval, parent.frame())
    if (any(sapply(args, length) > 1L))
        stop("all arguments must be logical(1)")
    if (length(args) > 0)
    {
        ## deprecate isValidVendorRead
        if ("isValidVendorRead" %in% names(args))
        {
            .Deprecated("isNotPassingQualityControls",
                        old="isValidVendorRead")
            old <- !args[["isValidVendorRead"]]
            args[["isValidVendorRead"]] <- NULL
            if (("isNotPassingQualityControls" %in% names(args)) &&
                !identical(args[["isNotPassingQualityControls"]], old))
            {
                msg <- sprintf("'%s' inconsistent with '%s', using '%s'",
                    "isValidVendorRead", "isNotPassingQualityControls",
                     "isNotPassingQualityControls")
                warning(paste(strwrap(msg, exdent=2), collapse="\n"))
            }
            args[["isNotPassingQualityControls"]] <- old
        }
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
    nms <- names(.scanBamTemplate())
    nms[nms != "tag"]
}

setMethod(show, "ScanBamParam",
          function(object)
{
    .show_classname(object)
    flag <- bamFlag(object)[!is.na(bamFlag(object))]
    vals <- sprintf("bamFlag (NA unless specified): %s",
                    paste(names(flag), flag, sep="=", collapse=", "))
    cat(strwrap(vals, exdent=2), sep="\n")
    cat("bamSimpleCigar: ", bamSimpleCigar(object), "\n", sep="")
    cat("bamReverseComplement: ", bamReverseComplement(object), "\n",
        sep="")
    cat("bamTag:", paste(bamTag(object), collapse=", "), "\n")
    cat("bamWhich:", length(bamWhich(object)), "elements\n")
    what <- paste("bamWhat: ", paste(bamWhat(object), collapse=", "))
    cat(strwrap(what, exdent=2), sep="\n")
})

## flag utils

## Explode the bits of a 'flag' vector into a matrix.
bamFlagAsBitMatrix <- function(flag, bitnames=FLAG_BITNAMES)
{
    bitpos <- match(bitnames, FLAG_BITNAMES)
    ans <- IRanges:::explodeIntBits(flag, bitpos=bitpos)
    dimnames(ans) <- list(names(flag), bitnames)
    ans
}

## Performs a logical AND between 2 'flag' vectors.
bamFlagAND <- function(flag1, flag2)
{
    bits1 <- bamFlagAsBitMatrix(flag1)
    bits2 <- bamFlagAsBitMatrix(flag2)
    ans <- IRanges:::implodeIntBits(bits1 & bits2)
    names(ans) <- names(flag1)
    ans
}

bamFlagTest <- function(flag, value)
{
    if (length(value) != 1 || !value %in% FLAG_BITNAMES) {
        msg <- sprintf("'is' must be character(1) in '%s'",
                       paste(FLAG_BITNAMES, collapse="' '"))
        stop(msg)
    }
    i <- 2 ^ (match(value, FLAG_BITNAMES) - 1L)
    bitAnd(flag, i) == i
}
