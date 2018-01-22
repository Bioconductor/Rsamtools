TabixFile <-
    function(file, index=paste(file, "tbi", sep="."), ...,
             yieldSize=NA_integer_)
{
    if (is(file, "TabixFile"))
        return(file)
    tryCatch({
        .io_check_exists(c(file, index))
    }, error=function(err) {
        stop(sprintf("TabixFile: %s", conditionMessage(err)), call.=FALSE)
    })
    .RsamtoolsFile(.TabixFile, file, index, yieldSize=yieldSize, ...)
}

open.TabixFile <-
    function(con, ...)
{
    ## FIXME: path? index?
    con$.extptr <- .Call(.tabixfile_open, path(con), index(con, asNA=FALSE))
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
        .Call(.index_tabix, file, format, as.integer(seq),
              as.integer(start), as.integer(end), skip, comment,
              zeroBased)
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

countTabix <-
    function(file, ...)
{
    tbxsym <- getNativeSymbolInfo(".tabix_count", "Rsamtools")
    scanTabix(file, ..., tbxsym=tbxsym)
}

.tabix_scan <-
    function(file, ..., space, start, end,
             tbxsym=getNativeSymbolInfo(".tabix_as_character",
               "Rsamtools"), tbxstate=NULL, row.names=NULL)
{
    tryCatch({
        if (!isOpen(file)) {
            open(file)
            on.exit(close(file))
        }

        if (0L != length(start) && !is.na(yieldSize(file))) {
            msg <- sprintf("'%s' must be '%s' when range specified",
                           "yieldSize(file)", "NA_integer_")
            cond <- simpleError(msg)
            class(cond) <- c("scanTabix_param", class(cond))
            stop(cond)
        }

        result <- .Call(.scan_tabix, .extptr(file),
                        list(space, start, end), yieldSize(file),
                        tbxsym$address, tbxstate, row.names)
        setNames(result, sprintf("%s:%d-%d", space, start, end))
    }, scanTabix_param=function(err) stop(err), error=function(err) {
        msg <- paste0("scanTabix: ", conditionMessage(err),
                      "\n path: ", path(file), "\n index: ",
                      index(file))
        cond <- simpleError(msg)
        class(cond) <- c("scanTabix_io", class(err))
        stop(cond)
    })
}

setMethod(scanTabix, c("TabixFile", "missing"),
    function(file, ..., param)
{
    .tabix_scan(file, ..., space=character(), start=integer(),
                end=integer())
})

setMethod(scanTabix, c("TabixFile", "IntegerRangesList"),
    function(file, ..., param)
{
    .tabix_scan(file, ..., space=as.character(space(param)),
                start=.uunlist(start(param)),
                end=.uunlist(end(param)))
})

setMethod(scanTabix, c("TabixFile", "GRanges"),
    function(file, ..., param)
{
    .tabix_scan(file, ..., space=as.character(seqnames(param)),
                start=start(param), end=end(param))
})

setMethod(scanTabix, c("character", "missing"),
    function(file, ..., param)
{
    scanTabix(TabixFile(file), ...)
})

setMethod(scanTabix, c("character", "ANY"),
    function(file, ..., param)
{
    scanTabix(TabixFile(file), ..., param=param)
})

setMethod(yieldTabix, "TabixFile",
    function(file, ..., yieldSize=1000000L)
{
    .Deprecated("scanTabix")
    scanTabix(file, ..., yieldSize=yieldSize)[[1]]
})
