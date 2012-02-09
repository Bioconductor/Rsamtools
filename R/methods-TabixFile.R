TabixFile <-
    function(file, index=paste(file, "tbi", sep="."), ...)
{
    tryCatch({
        .io_check_exists(c(file, index))
    }, error=function(err) {
        stop(sprintf("TabixFile: %s", conditionMessage(err)), call.=FALSE)
    })
    .RsamtoolsFile(.TabixFile, file, index)
}

open.TabixFile <-
    function(con, ...)
{
    ## FIXME: path? index?
    con$.extptr <- .Call(.tabixfile_open, path(con), index(con))
    invisible(con)
}

close.TabixFile <-
    function(con, ...)
{
    if (!isOpen(con))
        stop("isOpen(<TabixFile>) is not 'TRUE'")
    con$.extptr <- .Call(.tabixfile_close, .extptr(con))
    invisible(con)
}

setMethod(isOpen, "TabixFile",
    function(con, rw="")
{
    if (!missing(rw) && rw != "read")
        stop("'rw' must be 'read'")
    .Call(.tabixfile_isopen, .extptr(con))
})

indexTabix <- 
    function(file,
             format=c("gff", "bed", "sam", "vcf", "vcf4", "psltbl"),
             seq=integer(), start=integer(), end=integer(),
             skip=0L, comment="#", zeroBased=FALSE, ...)
{
    tryCatch({
        file <- .normalizePath(file)
        format <- 
            if (!missing(format)) match.arg(format)
            else character()
        idx <- .Call(.index_tabix, file, format,
                     as.integer(seq), as.integer(start),
                     as.integer(end), skip, comment, zeroBased)
        sprintf("%s.tbi", file)
    }, error=function(err) {
        stop(conditionMessage(err), "\n  file: ", file)
    })
}

.headerTabix <-
    function(file, ...)
{
    if (!isOpen(file)) {
        open(file)
        on.exit(close(file))
    }
    .Call(.header_tabix, .extptr(file))
}

setMethod(headerTabix, "TabixFile", .headerTabix)

setMethod(headerTabix, "character", function(file, ...) {
    .headerTabix(TabixFile(file))
})

.seqnamesTabix <-
    function(file, ...)
{
    .headerTabix(file, ...)[["seqnames"]]
}

setMethod(seqnamesTabix, "TabixFile", .seqnamesTabix)

setMethod(seqnamesTabix, "character", function(file, ...) {
    .seqnamesTabix(TabixFile(file))
})

.tabix_scan <-
    function(file, ..., space, start, end,
             tbxidx=rep(list(NULL), length(space)), yieldSize=30000L,
             tbxsym=getNativeSymbolInfo(".tabix_as_character",
               "Rsamtools"), tbxstate=NULL)
{
    tryCatch({
        if (!isOpen(file)) {
            open(file)
            on.exit(close(file))
        }

        tbxidx <- lapply(tbxidx, function(elt) {
            if (is.null(elt)) elt else sort(unique(as.integer(elt)))
        })
        if (length(tbxidx) != length(space))
            stop("length(tbxidx) must equal length(space(param))")
        result <- .Call(.scan_tabix, .extptr(file),
                        list(space, start, end), tbxidx,
                        as.integer(yieldSize), tbxsym$address,
                        tbxstate)
        names(result) <- sprintf("%s:%d-%d", space, start, end)
        result
    }, error=function(err) {
        stop("scanTabix: ", conditionMessage(err), "\n  path: ",
             path(file), call.=FALSE)
    })
}

setMethod(scanTabix, c("TabixFile", "RangesList"),
    function(file, ..., param)
{
    .tabix_scan(file, ..., space=as.character(space(param)),
                start=.uunlist(start(param)),
                end=.uunlist(end(param)))
})

setMethod(scanTabix, c("TabixFile", "RangedData"),
    function(file, ..., param)
{
    .tabix_scan(file, ..., param=ranges(param))
})

setMethod(scanTabix, c("TabixFile", "GRanges"),
    function(file, ..., param)
{
    .tabix_scan(file, ..., space=as.character(seqnames(param)),
                start=start(param), end=end(param))
})

setMethod(scanTabix, c("character", "RangesList"),
    function(file, ..., param)
{
    file <- TabixFile(file)
    callGeneric(file, ..., param=param)
})

setMethod(scanTabix, c("character", "RangedData"),
    function(file, ..., param)
{
    file <- TabixFile(file)
    callGeneric(file, ..., param=param)
})

setMethod(scanTabix, c("character", "GRanges"),
    function(file, ..., param)
{
    file <- TabixFile(file)
    callGeneric(file, ..., param=param)
})

.tabix_yield <-
    function(file, ..., yieldSize, tbxidx=NULL, tbxgrow=FALSE,
             tbxsym=getNativeSymbolInfo(".tabix_as_character",
               "Rsamtools"), tbxstate=NULL)

{
    tryCatch({
        if (!isOpen(file))
            open(file)
        if (!is.null(tbxidx))
            tbxidx <- sort(unique(as.integer(tbxidx)))
        .Call(.yield_tabix, .extptr(file), tbxidx,
              as.integer(yieldSize), as.logical(tbxgrow),
              tbxsym$address, tbxstate)
    }, error=function(err) {
        stop("yield: ", conditionMessage(err), "\n  path: ",
             path(file), call.=FALSE)
    })
}

setMethod(yieldTabix, "TabixFile",
    function(file, ..., yieldSize=1000000L)
{
    .tabix_yield(file, ..., yieldSize=yieldSize)
})

