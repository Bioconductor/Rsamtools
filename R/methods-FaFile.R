FaFile <-
    function(file, index=sprintf("%s.fai", file),
                   gzindex=sprintf("%s.gzi", file))
{
    ans <- .RsamtoolsFile(.FaFile, file, index)
    gzindex(ans) <- .normalizePath(gzindex)
    ans
}

setMethod(gzindex, "FaFile",
    function(object, asNA=TRUE)
{
    gzindex <- object$gzindex
    if (asNA && ((length(gzindex) == 0L) || !nzchar(gzindex)))
        NA_character_
    else
        gzindex
})

setReplaceMethod("gzindex", "FaFile",
    function(object, value)
{
    stopifnot(length(value) == 1L)
    object$gzindex <- as.character(value)
    object
})

setMethod(gzindex, "FaFileList",
    function(object, asNA=TRUE)
{
    sapply(object, gzindex, asNA=asNA)
})

setReplaceMethod("gzindex", "FaFileList",
    function(object, value)
{
    stopifnot(length(value) == length(path(object)))
    for (i in seq_along(object))
        gzindex(object[[i]]) <- value[i]
    object
})

setMethod(show, "FaFile", function(object) {
    cat("class:", class(object), "\n")
    cat(.ppath("path", path(object)))
    cat(.ppath("index", index(object)))
    cat(.ppath("gzindex", gzindex(object)))
    cat("isOpen:", isOpen(object), "\n")
    cat("yieldSize:", yieldSize(object), "\n")
})

open.FaFile <-
    function(con, ...)
{
    .io_check_exists(path(con))
    tryCatch({
        con$.extptr <- .Call(.fafile_open, path(con), index(con, asNA=FALSE),
                                                      gzindex(con, asNA=FALSE))
    }, error=function(err) {
        stop(conditionMessage(err), "\n  file: ", path(con))
    })
    invisible(con)
}

close.FaFile <-
    function(con, ...)
{
    if (!isOpen(con))
        stop("isOpen(<FaFile>) is not 'TRUE'")
    tryCatch({
        con$.extptr <- .Call(.fafile_close, .extptr(con))
    }, error=function(err) {
        stop(conditionMessage(err), "\n  file: ", path(con))
    })
    invisible(con)
}

setMethod(isOpen, "FaFile",
    function(con, rw="")
{
    if (!missing(rw) && rw != "read")
        stop("'rw' must be 'read'")
    tryCatch({
        .Call(.fafile_isopen, .extptr(con))
    }, error=function(err) {
        stop(conditionMessage(err), "\n  file: ", path(con))
    })
})

setMethod(indexFa, "FaFile",
    function(file, ...)
{
    tryCatch({
        file$index <- paste0(.Call(.index_fa, path(file)), ".fai")
    }, error=function(err) {
        stop(conditionMessage(err), "\n  file: ", path(file))
    })
    file
})

setMethod(scanFaIndex, "FaFile",
    function(file, ...)
{
    what <- list(character(), integer(), NULL, NULL, NULL)
    tbl <- scan(index(file, asNA=FALSE), what=what, quiet=TRUE)
    GRanges(tbl[[1]], IRanges(1, width=tbl[[2]]),
            seqlengths=structure(tbl[[2]], .Names=tbl[[1]]))
})

setMethod(scanFaIndex, "FaFileList",
    function(file, ..., as=c("GRangesList", "GRanges"))
{
    lst <- lapply(file, scanFaIndex, ...)
    switch(match.arg(as), GRanges={
        unique(unlist(GRangesList(lst), use.names=FALSE))
    }, GRangesList={
        GRangesList(lst)
    })
})

setMethod(countFa, "FaFile",
    function(file, ...)
{
    if (!isOpen(file)) {
        open(file)
        on.exit(close(file))
    }
    tryCatch({
        .Call(.n_fa, .extptr(file))
    }, error=function(err) {
        stop(conditionMessage(err), "\n  file: ", path(file))
    })
})

.scanFa <-
    function(file, param, ...,
        as=c("DNAStringSet", "RNAStringSet", "AAStringSet"))
{
    if (!isOpen(file)) {
        open(file)
        on.exit(close(file))
    }

    as <- match.arg(as)
    type <- sub("StringSet", "", as)
    base <- sub("Set", "", as)
    lkup <- Biostrings::get_seqtype_conversion_lookup("B", type)
    tryCatch({
        regions <- .asRegions(param)
        dna <- .Call(.scan_fa, .extptr(file), regions[[1]], regions[[2]],
                     regions[[3]], base, lkup)
        setNames(dna, regions[[1]])
    }, error=function(err) {
        stop(conditionMessage(err), "\n  file: ", path(file))
    })
}

setMethod(scanFa, c("FaFile", "GRanges"), .scanFa)

setMethod(scanFa, c("FaFile", "IntegerRangesList"), .scanFa)

setMethod(scanFa, c("FaFile", "missing"),
    function(file, param, ...,
        as=c("DNAStringSet", "RNAStringSet", "AAStringSet"))
{
    as <- match.arg(as)
    switch(as,
           DNAStringSet=readDNAStringSet(path(file), ...),
           RNAStringSet=readRNAStringSet(path(file), ...),
           AAStringSet=readAAStringSet(path(file), ...))
})

setMethod(seqinfo, "FaFile",
    function(x)
{
    gr <- scanFaIndex(x)
    Seqinfo(as.character(seqnames(gr)), width(gr))
})

setMethod(getSeq, "FaFile",
    function(x, param, ...)
{
    if (missing(param)) {
        scanFa(x, ...)
    } else {
        dna <- scanFa(x, param, ...)
        if (is(param, "GRanges")) {
            idx <- as.logical(strand(param) == "-")
            if (any(idx))
                dna[idx] <- reverseComplement(dna[idx])
        }
        dna
    }
})

setMethod(getSeq, "FaFileList",
    function(x, param, ...)
{
    if (!is(param, "GRangesList"))
        stop("'param' must be 'GRangesList' when 'x' is 'FaFileList'")
    if (length(x) != length(param))
        stop("'x' and 'param' must have equal length")

    lst <- lapply(seq_len(length(x)), function(i, x, param, ...) {
        getSeq(x[[i]], param[[i]], ...)
    }, x, param, ...)
    do.call(SimpleList, lst)
})

## character wrappers

setMethod(indexFa, "character",
    function(file, ...) index(indexFa(open(FaFile(file)))))

setMethod(scanFaIndex, "character",
    function(file, ...) scanFaIndex(open(FaFile(file))))

setMethod(countFa, "character",
    function(file, ...) countFa(FaFile(file)))

setMethod(scanFa, c("character", "GRanges"),
    function(file, param, ...,
        as=c("DNAStringSet", "RNAStringSet", "AAStringSet"))
{
    as <- match.arg(as)
    scanFa(FaFile(file), param, ..., as=as)
})

setMethod(scanFa, c("character", "IntegerRangesList"),
    function(file, param, ...,
        as=c("DNAStringSet", "RNAStringSet", "AAStringSet"))
{
    as <- match.arg(as)
    scanFa(FaFile(file), param, ..., as=as)
})

setMethod(scanFa, c("character", "missing"),
    function(file, param, ...,
        as=c("DNAStringSet", "RNAStringSet", "AAStringSet"))
{
    as <- match.arg(as)
    scanFa(FaFile(file), ..., as=as)
})
