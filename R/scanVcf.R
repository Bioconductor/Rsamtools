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
    function(fmt, tag, ...)
{
    ## map header FORMAT to R types
    map <- lapply(fmt$Type, switch,
                  String=character(), Integer=integer(),
                  Float=numeric(), Flag=logical())
    map[fmt$Number != 1] <- list(character())
    names(map) <- rownames(fmt)

    ## user selected 
    if (!identical(character(), tag))
        if (1L == length(tag) && is.na(tag)) {
            map[] <- list(NULL)
        } else {
            map[!names(map) %in% tag] <- list(NULL)
            if (!all(tag %in% names(map)))
                warning(paste("values in ScanVcfParam(", tag, "=c(...)) ",
                        "not present in header :", 
                        tag[!tag %in% names(map)], sep=""))
        }
    map
}

.vcf_scan <-
    function(file, ..., info=character(), geno=character(), 
             space, start, end)
{
    tryCatch({
        space <- as.character(space)    # could be factor
        if (!isOpen(file)) {
            open(file)
            on.exit(close(file))
        }
        hdr <- scanVcfHeader(file)[[1]]
        samples <- hdr$Sample
        imap <- .vcf_map(hdr$Header[["INFO"]], info)
        gmap <- .vcf_map(hdr$Header[["FORMAT"]], geno)
        yieldSize <- 1000000L           # a guess, grows as necessary
        result <- .Call(.scan_vcf, .extptr(file), list(space, start,
                        end), yieldSize, samples, imap, gmap)
        names(result) <- sprintf("%s:%d-%d", space, start, end)
        result
    }, error = function(err) {
        stop("scanVcf: ", conditionMessage(err), "\n  path: ", 
             path(file), call. = FALSE)
    })
}

.vcf_scan_connection <-
    function(file, ..., info=character(), geno=character())
{
    tryCatch({
        fl <- summary(file)$description
        hdr <- scanVcfHeader(fl)[[1]]
        samples <- hdr$Sample
        imap <- .vcf_map(hdr$Header[["INFO"]], info)
        gmap <- .vcf_map(hdr$Header[["FORMAT"]], geno)

        txt <- readLines(file, ...)
        txt <- txt[!grepl("^#", txt)] # FIXME: handle header lines better
        result <- .Call(.scan_vcf_connection, txt, samples, imap, gmap)
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
    ## no ranges
    if (length(vcfWhich(param)) == 0) 
        res <- callGeneric(path(file), param=param)
    else 
    ## ranges
        res <- scanVcf(file, ..., info=vcfInfo(param), geno=vcfGeno(param), 
                       param=vcfWhich(param))
    if (vcfTrimEmpty(param))
        lapply(res, function(rng) {
            rng[["GENO"]] <- Filter(Negate(is.null), rng[["GENO"]])
            rng
        })
    else res
})

setMethod(scanVcf, c("TabixFile", "missing"),
    function(file, ..., param)
{
    callGeneric(path(file))
})

setMethod(scanVcf, c("character", "ScanVcfParam"),
    function(file, ..., param)
{
    ## no ranges
    if (length(vcfWhich(param)) == 0) {
        con <- file(file)
        on.exit(close(con))
        .vcf_scan_connection(con, info=vcfInfo(param), geno=vcfGeno(param))
    } else {
    ## ranges
      callGeneric(TabixFile(file), ..., param=param)
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

.unpackVcfTag <-
    function(tag, id, n, type)
    ## 'tag' is a named list of info or geno,
    ## id, n, type are equal-length vectors of INFO or FORMAT info
{
    if (is.null(names(tag)))
        stop(print(tag, " must be a named list", sep=""))
    Map(function(elt, nm, id, n, type) {
        if (is.na(idx <- match(nm, id))) {
            msg <- sprintf("element '%s' not found in file header identifiers",
                           nm)
            stop(msg)
        }
    ## handle NULL elements
        if (is.null(elt))
            elt
        else 
            .unpackVcfField(elt, id[idx], n[idx], type[idx])
    }, tag, names(tag), MoreArgs=list(id, n, type))
}
 
setMethod(unpackVcf, c("list", "missing"),
    function(x, hdr, ..., info=TRUE, geno=TRUE)
{
    if (!is.logical(info))
        x <- lapply(x, function(elt, id, n, type) {
            elt[["INFO"]] <- .unpackVcfTag(elt[["INFO"]], id, n, type)
            elt
        }, rownames(info), info$Number, info$Type)
    if (!is.logical(geno))
        x <- lapply(x, function(elt, id, n, type) {
            elt[["GENO"]] <- .unpackVcfTag(elt[["GENO"]], id, n, type)
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
