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
        stop(conditionMessage(err), "\n  file: ", path(file))
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
    if (!missing(rw) && rw == "read")
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

setMethod(countFa, "FaFile",
    function(file, ...)
{
    tryCatch({
        .Call(.n_fa, .extptr(file))
    }, error=function(err) {
        stop(conditionMessage(err), "\n  file: ", path(file))
    })
})

setMethod(scanFa, c("FaFile", "GRanges"),
    function(file, param=GRanges(), ...)
{
    if (!isOpen(file))
        stop("'FaFile' not open")
    lkup <- Biostrings:::get_xsbasetypes_conversion_lookup("B", "DNA")
    tryCatch({
        nms <- as.character(seqnames(param))
        dna <- .Call(.scan_fa, .extptr(file), nms, start(param),
                     end(param), lkup)
        names(dna) <- nms
        dna
    }, error=function(err) {
        stop(conditionMessage(err), "\n  file: ", path(file))
    })
})

setMethod(scanFa, c("FaFile", "missing"),
    function(file, param=GRanges(), ...)
{
    if (!isOpen(file))
        stop("'FaFile' not open")
    read.DNAStringSet(path(file))
})

## character wrappers

setMethod(indexFa, "character",
    function(file, ...) index(indexFa(open(FaFile(file)))))

setMethod(scanFaIndex, "character",
    function(file, ...) scanFaIndex(open(FaFile(file))))

setMethod(countFa, "character",
    function(file, ...) countFa(open(FaFile(file))))

setMethod(scanFa, c("character", "GRanges"),
    function(file, param=GRanges(), ...)
{
    scanFa(open(FaFile(file)), param, ...)
})

setMethod(scanFa, c("character", "missing"),
    function(file, param, ...)
{
    scanFa(open(FaFile(file)), ...)
})
