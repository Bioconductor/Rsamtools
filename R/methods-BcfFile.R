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
    con$.extptr <-
        .Call(.bcffile_open, path(con), index(con, asNA=FALSE), bcfMode(con))
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
    h <- h[!duplicated(h)]

    ## Simple key=value pairs
    rex <- "^[[:alnum:]]+=[^<]"
    if (length(h1 <- grep(rex, h, value=TRUE))) {
        idx <- regexpr("=", h1, fixed=TRUE) ## first match
        rnms <- substring(h1, 1, idx - 1)
        if (is(rnms, "character")) {
            if (any(duplicated(rnms))) { 
                warning("duplicate keys in header will be forced to unique ",
                        "rownames")
                rnms <- make.unique(rnms)
            }
        }
        meta_df <- DataFrame(row.names=rnms, Value=substring(h1, idx + 1))
        meta <- as(splitAsList(meta_df, rnms), "SimpleDataFrameList")
    } else {
        meta <- DataFrameList()
    } 

    ## Non-simple key-value pairs 
    ## These lines have values enclosed in '<>' with subfields such as 'ID', 
    ## 'Number', 'Type' etc.
    rex <- "^([[:alnum:]]+)=<(.*)>"
    lines <- grep(rex, h, value=TRUE)
    if (length(lines) > 0) {
        if (any(duplicated(lines))) { 
            warning("duplicate INFO and FORMAT header lines will be ignored")
            lines <- lines[duplicated(lines) == FALSE]
        } 
    }
    tags <- sub(rex, "\\1", lines)

    keyval0 <- sub(rex, "\\2", lines)
    ## Handle INFO, FORMAT, FILTER, ALT, SAMPLE
    keyval1 <- rep(NA_character_, length(keyval0))
    keyval <- list()
    idx <- tags %in% c("INFO", "FORMAT", "FILTER", "ALT")
    keyval1[idx] <- strsplit(keyval0[idx], 
            ",(?=(ID|Number|Type)=[[:alnum:]]*)|,(?=Description=\".*?\")", 
            perl=TRUE)
    keyval[idx] <- lapply(which(idx), 
        function(i, keyval1) strsplit(keyval1[[i]],
            "(?<=[ID|Number|Type|Description])=", perl=TRUE), keyval1)
    keyval1[!idx] <- strsplit(keyval0[!idx], ",(?=[[:alnum:]]+=)", perl=TRUE)
    keyval[!idx] <- lapply(which(!idx), function(i, keyval1) {
        strsplit(keyval1[[i]], "(?<=[[:alnum:]])=", perl=TRUE)}, keyval1)

    tbls <- tapply(keyval, tags, 
        function(elt) {
            keys <- lapply(elt, sapply, "[[", 1)
            vals0 <- lapply(elt, sapply, "[[", 2)
            vals <- Map("names<-", vals0, keys)
            cols <- unique(unlist(keys))
            entries <- Map(function(k) as.vector(sapply(vals, "[", k)),
                           cols)
            desc <- which("DESCRIPTION" == toupper(names(entries)))
            if (1L == length(desc))
                entries[[desc]] <- gsub("\"", "", entries[[desc]])
            id <- which("ID" == toupper(names(entries)))
            if (length(id) > 0L) {
                if (any(duplicated(entries[[id]]))) 
                    warning("duplicate ID's in header will be forced to unique ",
                            "rownames")
                df <- DataFrame(entries[-id], 
                                row.names=make.unique(entries[[id]]))
            } else {
                ## ID is not a required field
                df <- DataFrame(entries)
            }
            df
        })

    ## 'GT' first in order
    if (length(tbls))
        if (length(tbls$FORMAT))
            if (any(GT <- rownames(tbls$FORMAT) %in% "GT"))
                tbls$FORMAT <- rbind(tbls$FORMAT[GT,], tbls$FORMAT[!GT,]) 

    SimpleList(Reference=header[["Reference"]],
               Sample=header[["Sample"]],
               Header=c(meta, tbls[unique(tags)]))
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
        geno <- bcfGeno(param)
        idx <- 
            if (1L == length(geno) && is.na(geno)) FALSE
            else names(tmpl[["GENO"]]) %in% geno
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
    regions <-
        if (0L != length(which)) {
            if (!nzchar(index(file)))
                stop("scanBcf 'index' must exist when non-zero bcfWhich()")
            list(as.character(space(which)), .uunlist(start(which)),
                 .uunlist(end(which)))
        } else NULL
    res <- tryCatch({
        .Call(func, .extptr(file), regions, tmpl)
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
    if (length(res[["GENO"]])) {
        idx <- sapply(res[["GENO"]], function(x) 2L == length(dim(x)))
        res[["GENO"]][idx] <- lapply(res[["GENO"]][idx], t)
    }
    res
}

setMethod(scanBcf, "BcfFile",
          function(file, ..., param=ScanBcfParam())
{
    .io_bcf(.scan_bcf, file, ..., param=param)
})

setMethod(asBcf, "character",
    function(file, dictionary, destination, ...,
             overwrite=FALSE, indexDestination=TRUE)
{
    file <- .normalizePath(file)
    destination <- .normalizePath(destination)
    destination <- paste(destination, "bcf", sep=".")

    dict <- tempfile()
    on.exit(unlink(dict))
    tryCatch({
        if (!overwrite && file.exists(destination)) {
            msg <-
                sprintf("'destination' exists, 'overwrite' is FALSE\n  destination.bcf: %s",
                        "destination", "overwrite", destination)
            stop(msg)
        }
        writeLines(dictionary, dict)
        result <- .Call(.as_bcf, file, dict, destination)
        if (!file.exists(destination))
            stop("failed to create 'BCF' file")
        if (indexDestination)
            indexBcf(destination)
    }, error=function(err) {
        msg <- sprintf("'asBcf' %s\n  VCF file: '%s'\n  destination: '%s'\n",
                       conditionMessage(err), file, destination)
        stop(msg)
    })
    destination
})

setMethod(indexBcf, "BcfFile", function(file, ...)
{
    if (!isOpen(file)) {
        open(file)
        on.exit(close(file))
    }
    if (.Call(.bcffile_isvcf, .extptr(file)))
        stop("'indexBcf' requires a BCF (not VCF) file")
    tryCatch({
        .Call(.index_bcf, path(file))
    }, error=function(err) {
        msg <- sprintf("'indexBcf' %s\n  file: '%s'\n",
                       conditionMessage(err), file)
        stop(msg)
    })
})

setMethod(indexBcf, "character",
          function(file, ...)
{
    indexBcf(BcfFile(file, character()), ...)
})
