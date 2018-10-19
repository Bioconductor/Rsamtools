## Native method first.
setMethod(ScanBamParam, c(which="IntegerRangesList"),
          function(flag=scanBamFlag(), simpleCigar=FALSE,
                   reverseComplement=FALSE, tag=character(0),
                   tagFilter=list(), what=character(0), which,
                   mapqFilter=NA_integer_)
{
    if (is.null(names(which))) {
        if (length(which) != 0L)
            stop(wmsg("'which' must have names in \"ScanBamParam\" ",
                      "method for IntegerRangesList objects "))
        names(which) <- character()
    }
    new("ScanBamParam", flag=flag, simpleCigar=simpleCigar,
        reverseComplement=reverseComplement, tag=tag,
        tagFilter=.normalize_tagFilter(tagFilter), what=what,
        which=which,
        mapqFilter=as.integer(mapqFilter))
})

setMethod(ScanBamParam, c(which="missing"),
          function(flag=scanBamFlag(), simpleCigar=FALSE,
                   reverseComplement=FALSE, tag=character(0),
                   tagFilter=list(), what=character(0), which,
                   mapqFilter=NA_integer_)
{
    which <- IRangesList()
    ScanBamParam(flag=flag, simpleCigar=simpleCigar,
                 reverseComplement=reverseComplement, tag=tag,
                 tagFilter=tagFilter, what=what, which=which,
                 mapqFilter=as.integer(mapqFilter))
})

## Default method.
setMethod(ScanBamParam, c(which="ANY"),
          function(flag=scanBamFlag(), simpleCigar=FALSE,
                   reverseComplement=FALSE, tag=character(0),
                   tagFilter=list(), what=character(0), which,
                   mapqFilter=NA_integer_)
{
    which <- as(which, "IntegerRangesList")
    ScanBamParam(flag=flag, simpleCigar=simpleCigar,
                 reverseComplement=reverseComplement, tag=tag,
                 tagFilter=tagFilter, what=what, which=which,
                 mapqFilter=as.integer(mapqFilter))
})

## Note that this method is not needed. Coercing a GRanges to IntegerRangesList
## works and produces exactly what the method below is doing by hand. So the
## default ScanBamParam method above just works on a GRanges. -- Herv\'e
setMethod(ScanBamParam, c(which="GRanges"),
          function(flag=scanBamFlag(), simpleCigar=FALSE,
                   reverseComplement=FALSE, tag=character(0),
                   tagFilter=list(), what=character(0), which,
                   mapqFilter=NA_integer_)
{
    which <- split(ranges(which), seqnames(which))
    ScanBamParam(flag=flag, simpleCigar=simpleCigar,
                 reverseComplement=reverseComplement, tag=tag,
                 tagFilter=tagFilter, what=what, which=which,
                 mapqFilter=as.integer(mapqFilter))
})

## adapted from ?integer
.has.wholenumbers <- function(x, tol=.Machine$double.eps^0.5)
    all(abs(x - round(x)) < tol)

.normalize_tagFilter <- function(tagfilter) {
    tagfilter <- lapply(tagfilter, function(x) {
        if(is.numeric(x)) {
            if(!.has.wholenumbers(x)) {
                msg <- paste0("Filtering tags by floating point values ",
                              "not supported. Please use integer values.")
                stop(paste(strwrap(msg, exdent=2), collapse="\n"))
            }
            x <- as.integer(x)
        }
        x ##unique(x)
    })
    tagfilter
}

setValidity("ScanBamParam", function(object) {
    msg <- NULL
    flag <- bamFlag(object, asInteger=TRUE)
    simpleCigar <- bamSimpleCigar(object)
    reverseComplement <- bamReverseComplement(object)
    tag <- bamTag(object)
    tagFilter <- bamTagFilter(object)
    what <- bamWhat(object)
    which <- bamWhich(object)
    mapqFilter <- bamMapqFilter(object)
    ## flag
    if (length(flag) != 2 || typeof(flag) != "integer")
        msg <- c(msg, "'flag' must be integer(2)")
    else {
        if (!identical(names(flag), c("keep0", "keep1")))
            msg <- c(msg, "'names(flag)' must be c('keep0', 'keep1')")
        else if (any(flag < 0 | flag >= 2^12))
            msg <- c(msg, "'flag' values must be >=0, <4096")
    }
    ## simpleCigar
    if (!((1L == length(simpleCigar)) &&
          !is.na(simpleCigar)))
        msg <- c(msg, "'simpleCigar' must be logical(1), not NA")
    ## reverseComplement
    if (!((1L == length(reverseComplement)) &&
          !is.na(reverseComplement)))
        msg <- c(msg, "'reverseComplement' must be logical(1)")
    ## tag
    if (0 != length(tag) && any(2 != nchar(tag)))
        msg <- c(msg, "'tag' must be two letters, e.g,, 'MD'")
    ## tagFilter
    if (length(tagFilter) > 0) {
        tagnms <- names(tagFilter)
        if(length(tagnms) == 0L || any(2 != nchar(tagnms)))
            msg <- c(msg, paste0("'tagFilter' elements must be named with ",
                                 "two characters each"))
        for(valvec in tagFilter)
            if(is.character(valvec) && any(!nzchar(valvec)))
                msg <- c(msg, "character() tagFilter values must be non-empty")
        elt_typeinfo <- vapply(tagFilter, function(x) {
            length(x) && (is.character(x) || is.integer(x)) && ##is.atomic(x)
            !is.null(x) && !anyNA(x)
        }, logical(1))
        if(any(!elt_typeinfo))
            msg <- c(msg, paste0("'tagFilter' must contain only non-NULL, ",
                                 "non-NA, non-empty character or integer ",
                                 "values"))
    }
    ## what
    if (!all(what %in% scanBamWhat()))
        msg <- c(msg, "'what' must be from 'scanBamWhat()'")
    ## which
    if (length(which) != 0 &&
        (any(!nzchar(names(which))) || any(is.na(names(which)))))
        msg <- c(msg, "'which' elements must be named (not NA)")
    ## mapqFilter
    if (length(mapqFilter) != 1)
        msg <- c(msg, "'mapqFilter' must be integer(1), >= 0")
    else if (!(is.na(mapqFilter) || mapqFilter >= 0))
        msg <- c(msg, "'mapqFilter' must be NA or >= 0")

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
    validObject(object)
    object
}
bamSimpleCigar <- function(object) slot(object, "simpleCigar")
"bamSimpleCigar<-" <- function(object, value)
{
    slot(object, "simpleCigar") <- value
    validObject(object)
    object
}
bamReverseComplement <-
    function(object) slot(object, "reverseComplement")
"bamReverseComplement<-" <- function(object, value)
{
    slot(object, "reverseComplement") <- value
    validObject(object)
    object
}
bamTag <- function(object) slot(object, "tag")
"bamTag<-" <- function(object, value)
{
    slot(object, "tag") <- value
    validObject(object)
    object
}
bamTagFilter <- function(object) slot(object, "tagFilter")
"bamTagFilter<-" <- function(object, value)
{
    slot(object, "tagFilter") <- .normalize_tagFilter(value)
    validObject(object)
    object
}

bamWhich <- function(object) slot(object, "which")

setGeneric("bamWhich<-",
           function(object, value) standardGeneric("bamWhich<-"))
        
setReplaceMethod("bamWhich", c("ScanBamParam", "IntegerRangesList"),
    function(object, value) 
{
    slot(object, "which") <- value
    validObject(object)
    object
})

setReplaceMethod("bamWhich", c("ScanBamParam", "GRanges"),
    function(object, value) 
{
    bamWhich(object) <- split(ranges(value), seqnames(value))
    validObject(object)
    object
})
        
setReplaceMethod("bamWhich", c("ScanBamParam", "ANY"),
    function(object, value) 
{
    bamWhich(object) <-  as(value, "IntegerRangesList")
    validObject(object)
    object
})

bamWhat <- function(object) slot(object, "what")
"bamWhat<-" <- function(object, value)
{
    slot(object, "what") <- value
    validObject(object)
    object
}

bamMapqFilter <- function(object) slot(object, "mapqFilter")
"bamMapqFilter<-" <- function(object, value)
{
    slot(object, "mapqFilter") <- as.integer(value)
    validObject(object)
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
    "isSecondaryAlignment",
    "isNotPassingQualityControls",
    "isDuplicate",
    "isSupplementaryAlignment"
)

scanBamFlag <-
    function(isPaired=NA, isProperPair=NA, isUnmappedQuery=NA,
             hasUnmappedMate=NA, isMinusStrand=NA,
             isMateMinusStrand=NA, isFirstMateRead=NA,
             isSecondMateRead=NA, # redundant
             isNotPrimaryRead=NA, isSecondaryAlignment=NA,
             isNotPassingQualityControls=NA, isDuplicate=NA,
             isSupplementaryAlignment = NA)

    ## NA: keep either 0 or 1 flag; FALSE: keep 0 flag; TRUE: keep 1 flag
{
    flag <- S4Vectors:::makePowersOfTwo(length(FLAG_BITNAMES))
    names(flag) <- FLAG_BITNAMES
    args <- lapply(as.list(match.call())[-1], eval, parent.frame())
    if (any(sapply(args, length) > 1L))
        stop("all arguments must be logical(1)")
    if (length(args) > 0)
    {
        ## deprecate isNotPrimaryRead
        if ("isNotPrimaryRead" %in% names(args))
        {
            .Deprecated("isSecondaryAlignment",
                        old="isNotPrimaryRead")
            old <- args[["isNotPrimaryRead"]]
            args[["isNotPrimaryRead"]] <- NULL
            value_to_use <- if ("isSecondaryAlignment" %in% names(args))
                args[["isSecondaryAlignment"]]
            else
                old
            if (("isSecondaryAlignment" %in% names(args)) &&
                !identical(args[["isSecondaryAlignment"]], old))
            {
                msg <- sprintf("'%s' inconsistent with '%s', using '%s'",
                               "isNotPrimaryRead", "isSecondaryAlignment",
                               "isSecondaryAlignment")
                warning(paste(strwrap(msg, exdent=2), collapse="\n"))
            }
            args[["isSecondaryAlignment"]] <- value_to_use
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
    cat("bamTagFilter:\n")
    tagFilter <- lapply(bamTagFilter(object), function(x) {
        if(is.character(x)) sQuote(x) else x
    })
    tagnms <- names(tagFilter)
    for(idx in seq_along(tagFilter)) {
        elt <- paste(tagnms[[idx]],
                     paste(S4Vectors:::selectSome(tagFilter[[idx]]),
                           collapse=", "), sep=" : ")
        cat(strwrap(elt, exdent=7, indent=2), sep="\n")
    }
    cat("bamWhich:", sum(lengths(bamWhich(object))), "ranges\n")
    what <- paste("bamWhat: ", paste(bamWhat(object), collapse=", "))
    cat(strwrap(what, exdent=2), sep="\n")
    cat("bamMapqFilter: ", bamMapqFilter(object), "\n", sep="")
})

## flag utils

## Explode the bits of a 'flag' vector into a matrix.
bamFlagAsBitMatrix <- function(flag, bitnames=FLAG_BITNAMES)
{
    ## deprecate 'isNotPrimaryRead'
    oldname <- "isNotPrimaryRead"
    newname <- "isSecondaryAlignment"
    if(oldname %in% bitnames && newname %in% bitnames) {
        msg <- paste0("'%s' is deprecated form of '%s', please specify ",
                      "only one")
        msg <- sprintf(msg, oldname, newname)
        stop(paste(strwrap(msg, exdent=2), collapse="\n"))
    }
    if(oldname %in% bitnames) {
        .Deprecated(newname, old=oldname)
        normalized <- bitnames
        normalized[which(normalized == oldname)] <- newname
        bitpos <- .calcBitPos(normalized)
    } else {
        bitpos <- .calcBitPos(bitnames)
    }
    ans <- S4Vectors:::explodeIntBits(flag, bitpos=bitpos)
    dimnames(ans) <- list(names(flag), bitnames)
    ans
}

.calcBitPos <- function(bitnames) {
    bitpos <- match(bitnames, FLAG_BITNAMES)
    invalid_bitnames_idx <- which(is.na(bitpos))
    if (length(invalid_bitnames_idx) != 0L) {
        in1string <- paste0(bitnames[invalid_bitnames_idx], collapse=", ")
        stop("invalid bitname(s): ", in1string)
    }
    bitpos
}

## Performs a logical AND between 2 'flag' vectors.
bamFlagAND <- function(flag1, flag2)
{
    bits1 <- bamFlagAsBitMatrix(flag1)
    bits2 <- bamFlagAsBitMatrix(flag2)
    ans <- S4Vectors:::implodeIntBits(bits1 & bits2)
    names(ans) <- names(flag1)
    ans
}

bamFlagTest <- function(flag, value)
{
    ## deprecate 'isNotPrimaryRead'
    oldname <- "isNotPrimaryRead"
    newname <- "isSecondaryAlignment"
    if(oldname == value) {
        .Deprecated(newname, old=oldname)
        value <- newname
    }
    if (length(value) != 1 || !value %in% FLAG_BITNAMES) {
        msg <- sprintf("'is' must be character(1) in '%s'",
                       paste(FLAG_BITNAMES, collapse="' '"))
        stop(msg)
    }
    i <- 2 ^ (match(value, FLAG_BITNAMES) - 1L)
    bitAnd(flag, i) == i
}
