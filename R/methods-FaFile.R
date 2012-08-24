FaFile <-
    function(file, ...)
{
    .RsamtoolsFile(.FaFile, file, file)
}

open.FaFile <-
    function(con, ...)
{
    .io_check_exists(path(con))
    tryCatch({
        con$.extptr <- .Call(.fafile_open, path(con))
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
        file$index <- .Call(.index_fa, path(file))
    }, error=function(err) {
        stop(conditionMessage(err), "\n  file: ", path(file))
    })
    file
})

setMethod(scanFaIndex, "FaFile",
    function(file, ...)
{
    what <- list(character(), integer(), NULL, NULL, NULL)
    tbl <- scan(sprintf("%s.fai", index(file)), what=what, quiet=TRUE)
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
    tryCatch({
        .Call(.n_fa, .extptr(file))
    }, error=function(err) {
        stop(conditionMessage(err), "\n  file: ", path(file))
    })
})

.scanFa <-
    function(file, param, ...)
{
    if (!isOpen(file)) {
        open(file)
        on.exit(close(file))
    }

    lkup <- Biostrings:::get_seqtype_conversion_lookup("B", "DNA")
    tryCatch({
        spc <- .asSpace(param)
        dna <- .Call(.scan_fa, .extptr(file), spc[[1]], spc[[2]],
                     spc[[3]], lkup)
        setNames(dna, spc[[1]])
    }, error=function(err) {
        stop(conditionMessage(err), "\n  file: ", path(file))
    })
}

setMethod(scanFa, c("FaFile", "GRanges"), .scanFa)

setMethod(scanFa, c("FaFile", "RangesList"), .scanFa)

setMethod(scanFa, c("FaFile", "RangedData"), .scanFa)

setMethod(scanFa, c("FaFile", "missing"),
    function(file, param, ...)
{
    readDNAStringSet(path(file))
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
    function(file, ...) countFa(open(FaFile(file))))

setMethod(scanFa, c("character", "GRanges"),
    function(file, param, ...)
{
    scanFa(FaFile(file), param, ...)
})

setMethod(scanFa, c("character", "RangesList"),
    function(file, param, ...)
{
    scanFa(FaFile(file), param, ...)
})

setMethod(scanFa, c("character", "RangedData"),
    function(file, param, ...)
{
    scanFa(FaFile(file), param, ...)
})

setMethod(scanFa, c("character", "missing"),
    function(file, param, ...)
{
    scanFa(FaFile(file), ...)
})
