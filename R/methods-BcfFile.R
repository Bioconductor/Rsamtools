BcfFile <-
    function(file, index=file,
             mode=ifelse(grepl("\\.bcf$", file), "rb", "r"))
{
    bf <- .RsamtoolsFile(.BcfFile, file, index)
    bf$mode <- mode
    bf
}

bcfMode <- function(object) object$mode

open.BcfFile <-
    function(con, ...)
{
    .io_check_exists(path(con))
    con$.extptr <- .Call(.bcffile_open, path(con), index(con), bcfMode(con))
    invisible(con)
}

close.BcfFile <-
    function(con, ...)
{
    if (!isOpen(con))
        stop("isOpen(<BcfFile>) is not 'TRUE'")
    con$.extptr <- .Call(.bcffile_close, .extptr(con))
    invisible(con)
}

setMethod(isOpen, "BcfFile",
    function(con, rw="")
{
    if (!missing(rw) && rw != "read")
        stop("'rw' must be 'read'")
    .Call(.bcffile_isopen, .extptr(con))
})

## scanBcfHeader

.bcfHeaderAsSimpleList <-
    function(header)
{
    h <- sub("^##", "", header[["Header"]])
    ## simple key=value pairs --> DataFrame
    rex <- "^[[:alpha:]]+=[^<]"
    x <- strsplit(grep(rex, h, value=TRUE), "=")
    meta <- DataFrame(row.names=sapply(x, "[[", 1),
                      Value=sapply(x, "[[", 2))
    ## key=<values> as SimpleList of DataFrame's
    rex <- "^([[:alpha:]]+)=<(.*)>"
    lines <- grep(rex, h, value=TRUE)
    tags <- sub(rex, "\\1", lines)
    vals <- strsplit(sub(rex, "\\2", lines), "[=,]")
    tbls <- tapply(vals, tags, function(elt) {
        keys <- lapply(elt, "[", c(TRUE, FALSE))
        vals0 <- lapply(elt, "[", c(FALSE, TRUE))
        vals <- Map("names<-", vals0, keys)
        cols <- unique(unlist(keys))
        entries <- Map(function(k) as.vector(sapply(vals, "[", k)),
                       cols)
        desc <- grep("DESCRIPTION", toupper(names(entries)))
        if (1L == length(desc))
            entries[[desc]] <- gsub("\"", "", entries[[desc]])
        id <- grep("ID", toupper(names(entries)))
        df <- do.call(DataFrame, entries[-id])
        rownames(df) <- entries[[id]]
        df
    })
    SimpleList(Reference=header[["Reference"]],
               Sample=header[["Sample"]],
               Header=c(DataFrameList(META=meta), tbls))
}

setMethod(scanBcfHeader, "BcfFile",
          function(file, ...)
{
    if (!isOpen(file)) {
        open(file)
        on.exit(close(file))
    } else if ("r" == bcfMode(file)) {
        on.exit(open(file))             # re-open
    }
    header <- .Call(.scan_bcf_header, .extptr(file))
    .bcfHeaderAsSimpleList(header)
})

## scanBcf

.bcf_template <-
    function(param)
{
    tmpl <- list(CHROM=character(),
                 POS=integer(),
                 ID=character(),
                 REF=character(),
                 ALT=character(),
                 QUAL=numeric(),
                 FILTER=character(),
                 INFO=character(), 
                 FORMAT=character(),
                 GENO=list(
                   PL=list(),
                   DP=structure(integer(0), .Dim=c(0L, 0L)),
                   GQ=structure(integer(0), .Dim=c(0L, 0L)),
                   SP=structure(integer(0), .Dim=c(0L, 0L)),
                   GT=structure(character(0), .Dim=c(0L, 0L)),
                   GL=list()),
                 RecordsPerRange=integer())
    if (!identical(character(), bcfGeno(param))) {
        idx <- names(tmpl[["GENO"]]) %in% bcfGeno(param)
        tmpl[["GENO"]] <- tmpl[["GENO"]][idx]
    }
    tmpl
}

.io_bcf <-
    function(func, file, ..., param)
{
    if (!isOpen(file)) {
        open(file)
        on.exit(close(file))
    } else if ("r" == bcfMode(file)) {
        on.exit(open(file))             # re-open
    }
    tmpl <- .bcf_template(param)
    which <- bcfWhich(param)
    space <- 
        if (0L != length(which)) {
            if (!nzchar(index(file)))
                stop("scanBcf 'index' must exist when non-zero bcfWhich()")
            list(as.character(space(which)), .uunlist(start(which)),
                 .uunlist(end(which)))
        } else NULL
    res <- tryCatch({
        .Call(func, .extptr(file), space, tmpl)
    }, error=function(err) {
        stop("scanBcf: ", conditionMessage(err), "\n  path: ",
             path(file), call.=FALSE)
    })
    if (bcfTrimEmpty(param)) {
        idx <- sapply(res[["GENO"]], function(x) {
            tst <- switch(typeof(x), list=is.null, is.na)
            !all(unlist(lapply(x, tst), use.names=FALSE))
        })
        res[["GENO"]] <- res[["GENO"]][idx]
    }
    idx <- sapply(res[["GENO"]], function(x) 2L == length(dim(x)))
    res[["GENO"]][idx] <- lapply(res[["GENO"]][idx], t)
    res
}

setMethod(scanBcf, "BcfFile",

          function(file, ..., param=ScanBcfParam())
{
    res <- .io_bcf(.scan_bcf, file, ..., param=param)
})
