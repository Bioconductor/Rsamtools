## scanVcfHeader

setMethod(scanVcfHeader, "character",
    function(file, ...) 
{
    scanBcfHeader(file, ...)
})

setMethod(scanVcfHeader, "TabixFile",
    function(file, ...)
{
    scanVcfHeader(path(file), ...)
})


## scanVcf 

.vcf_map <-
    function(hdr, geno, ...)
{
    ## map header FORMAT to R types
    fmt <- hdr$Header[["FORMAT"]]
    map <- lapply(fmt$Type, switch,
                  String=character(), Integer=integer(),
                  Float=numeric())
    map[fmt$Number != 1] <- list(character())
    names(map) <- rownames(fmt)

    ## user selected geno
    if (!identical(character(), geno))
        if (1L == length(geno) && is.na(geno)) {
            map[] <- list(NULL)
        } else {
            map[!names(map) %in% geno] <- list(NULL)
            if (!all(geno %in% names(map)))
                warning("values in ScanVcfParam(geno=c(...)) ",
                        "not present in FORMAT field :", 
                        geno[!geno %in% names(map)])
        }
    map
}

.vcf_scan <-
    function(file, ..., geno=character(), space, start, end)
{
    tryCatch({
        space <- as.character(space)    # could be factor
        if (!isOpen(file)) {
            open(file)
            on.exit(close(file))
        }
        hdr <- scanVcfHeader(file)[[1]]
        samples <- hdr$Sample
        map <- .vcf_map(hdr, geno=geno)
        yieldSize <- 1000000L           # a guess, grows as necessary
        result <- .Call(.scan_vcf, .extptr(file), list(space, start,
                        end), yieldSize, samples, map)
        names(result) <- sprintf("%s:%d-%d", space, start, end)
        result
    }, error = function(err) {
        stop("scanVcf: ", conditionMessage(err), "\n  path: ", 
             path(file), call. = FALSE)
    })
}

.vcf_scan_connection <-
    function(file, ..., geno=character())
{
    tryCatch({
        fl <- summary(file)$description
        hdr <- scanVcfHeader(fl)[[1]]
        samples <- hdr$Sample
        map <- .vcf_map(hdr, geno=geno)

        txt <- readLines(file, ...)
        txt <- txt[!grepl("^#", txt)] # FIXME: handle header lines better
        result <- .Call(.scan_vcf_connection, txt, samples, map)
        names(result) <- "*:*-*"
        result
    }, error = function(err) {
        stop("scanVcf: ", conditionMessage(err), "\n  path: ", 
             summary(file)$description, call. = FALSE)
    })
}

setMethod(scanVcf, c("TabixFile", "RangesList"),
    function(file, ..., param)
{
    .vcf_scan(file, ..., space=space(param),
              start=.uunlist(start(param)), end=.uunlist(end(param)))
})

setMethod(scanVcf, c("TabixFile", "RangedData"),
    function(file, ..., param)
{
    .vcf_scan(file, ..., param=ranges(param))
})

setMethod(scanVcf, c("TabixFile", "GRanges"),
    function(file, ..., param)
{
    .vcf_scan(file, ..., space=as.character(seqnames(param)),
              start=start(param), end=end(param))
})

setMethod(scanVcf, c("TabixFile", "ScanVcfParam"),
    function(file, ..., param)
{
    res <- scanVcf(file, ..., geno=vcfGeno(param), param=vcfWhich(param))
    if (vcfTrimEmpty(param))
        lapply(res, function(rng) {
            rng[["GENO"]] <- Filter(Negate(is.null), rng[["GENO"]])
            rng
        })
    else res
})

setMethod(scanVcf, c("character", "ANY"),
    function(file, ..., param)
{
    ## ScanVcfParam with no ranges
    if (class(param) == "ScanVcfParam") {
        if (length(vcfWhich(param)) == 0) {
            con <- file(file)
            on.exit(close(con))
            .vcf_scan_connection(con, geno=vcfGeno(param))
        } else {
    ## ScanVcfParam with ranges
          file <- TabixFile(file)
          scanVcf(file, ..., param=param)
        }
    } else {
    ## all others 
        file <- TabixFile(file)
        scanVcf(file, ..., param=param)
    }
})

setMethod(scanVcf, c("character", "missing"),
    function(file, ..., param)
{
    con <- file(file)
    on.exit(close(con))
    callGeneric(con, ...)
})

setMethod(scanVcf, c("connection", "missing"),
    function(file, ..., param)
{
    .vcf_scan_connection(file, ...)
})

## unpackVcf

.unpackVcfField <-
    function(x, id, n, type)
    ## convert sub-fields in 'x' to full R representation
{
    d <- suppressWarnings(as.integer(n))
    na <- all(is.na(unlist(x, use.names=FALSE)))
    withCallingHandlers({
        x <- if (na == TRUE) {
            x <- unlist(x, use.names=FALSE)
        } else {
            x <- if (!is.na(d)) {
                if (1L < d) {
                    idx <- as.integer(!is.na(x))
                    idx[idx == 0] <- d
                    xrep <- rep(x, idx)
                    x <- array(unlist(strsplit(xrep, ",", fixed=TRUE)),
                               dim=c(d, nrow(x), ncol(x)),
                               dimnames=list(NULL, NULL, colnames(x)))
                    x <- aperm(x, c(2, 3, 1))
                }
                switch(type, Flag= !is.na(x),
                       Character=, String=x,
                       Integer={ mode(x) <- "integer"; x },
                       Float={ mode(x) <- "numeric"; x },
                       stop(sprintf("unhandled FORMAT type '%s'", type)))
            } else {
                x <- split(strsplit(x, ",", fixed=TRUE), seq_len(nrow(x)))
                x <- switch(type, Character=, String=x,
                            Integer=lapply(x, lapply, as.integer),
                            Float=lapply(x, lapply, as.numeric),
                            stop(sprintf("unhandled FORMAT type '%s'", type)))
                names(x) <- NULL
                x
            }
        }
    }, warning=function(w) {
        msg <- sprintf("unpackVcf field '%s': %s", id,
                       conditionMessage(w))
        warning(msg, call.=FALSE)
        invokeRestart("muffleWarning")
    })
}

.unpackVcfInfo <-
    function(info, id, n, type)
{
    nrec <- length(info)

    recs <- strsplit(info, ";", fixed=TRUE)
    ridx <- rep(seq_along(recs), sapply(recs, length))
    flds <- strsplit(unlist(recs), "=", fixed=TRUE)
    tags <- factor(sapply(flds, "[[", 1), levels=id)
    names(type) <- id

    info <- Map(function(type, ridx, flds) {
        res <- matrix(NA_character_, nrow=nrec, ncol=1)
        res[ridx] <- 
            if (type=="Flag") ""
            else sapply(flds, "[[", 2)
        res
    }, type[levels(tags)], split(ridx, tags), split(flds, tags))

    result <- Map(.unpackVcfField, info, id, n, type)
    lapply(result, function(elt) {
        if (is(elt, "list"))
            unlist(elt, recursive=FALSE, use.names=FALSE)
        else elt
    })
}

.unpackVcfGeno <-
    function(geno, id, n, type)
    ## 'geno' is a named list,
    ## id, n, type are equal-length vectors of FORMAT info
{
    if (is.null(names(geno)))
        stop("'GENO' must be a named list")
    Map(function(elt, nm, id, n, type) {
        if (is.na(idx <- match(nm, id))) {
            msg <- sprintf("element '%s' not found in FORMAT identifiers",
                           nm)
            stop(msg)
        }
    ## handle NULL elements
        if (is.null(elt))
            elt
        else 
            .unpackVcfField(elt, id[idx], n[idx], type[idx])
    }, geno, names(geno), MoreArgs=list(id, n, type))
}
 
setMethod(unpackVcf, c("list", "missing"),
    function(x, hdr, ..., info=TRUE, geno=TRUE)
{
    if (!is.logical(info))
        x <- lapply(x, function(elt, id, n, type) {
            elt[["INFO"]] <- .unpackVcfInfo(elt[["INFO"]], id, n, type)
            elt
        }, rownames(info), info$Number, info$Type)
    if (!is.logical(geno))
        x <- lapply(x, function(elt, id, n, type) {
            elt[["GENO"]] <- .unpackVcfGeno(elt[["GENO"]], id, n, type)
            elt
        }, rownames(geno), geno$Number, geno$Type)
    x
})

setMethod(unpackVcf, c("list", "character") ,
    function(x, hdr, ..., info=TRUE, geno=TRUE)
{
    if (info) {
        if (length(x[[1]]$INFO) != 0) {
            info <- scanVcfHeader(hdr)[[1]][["Header"]][["INFO"]]
            if (is.null(info)) {
                msg <- sprintf("'INFO' vcf header info not found\n  path: %s",
                               hdr)
                warning(msg)
            info <- FALSE
            }
        } else {
            info <- FALSE
        }
    }
    if (geno) {
        if (length(unlist(x[[1]]$GENO, use.names=FALSE)) != 0) {
            geno <- scanVcfHeader(hdr)[[1]][["Header"]][["FORMAT"]]
            if (is.null(geno)) {
                msg <- sprintf("'FORMAT' vcf header info not found\n  path: %s",
                               hdr)
                warning(msg)
            geno <- FALSE
            }
       } else {
            geno <- FALSE
       }
    }
    unpackVcf(x, ..., info=info, geno=geno)
})

setMethod(unpackVcf, c("list", "TabixFile") ,
    function(x, hdr, ..., info=TRUE, geno=TRUE)
{
    unpackVcf(x, path(hdr), ..., info=info, geno=geno)
})
