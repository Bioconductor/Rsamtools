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
    function(hdr)
{
    ## map header FORMAT to R types
    fmt <- hdr$Header[["FORMAT"]]
    map <- lapply(fmt$Type, switch,
                  String=character(), Integer=integer(),
                  Float=numeric())
    map[fmt$Number != 1] <- list(character())
    names(map) <- rownames(fmt)
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
        map <- .vcf_map(hdr)
        if (!identical(character(), geno))
            if (1L == length(geno) && is.na(geno))
                map[] <- list(NULL)
            else
                map[!names(map) %in% geno] <- list(NULL)

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
    res <- scanVcf(file, ..., info=vcfInfo(param),
                   geno=vcfGeno(param), param=vcfWhich(param))
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
    file <- TabixFile(file)
    scanVcf(file, ..., param=param)
})

.vcf_scan_connection <-
    function(file, ...)
{
    tryCatch({
        fl <- summary(file)$description
        hdr <- scanVcfHeader(fl)[[1]]
        samples <- hdr$Sample
        map <- .vcf_map(hdr)

        txt <- readLines(file, ...)
        txt <- txt[!grepl("^#", txt)] # FIXME: handle header lines better
        result <- .Call(.scan_vcf_connection, txt, samples, map)
        names(result) <- fl
        result
    }, error = function(err) {
        stop("scanVcf: ", conditionMessage(err), "\n  path: ", 
             summary(file)$description, call. = FALSE)
    })
}


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

## utilities

.unpackVcfField <-
    function(x, id, n, type)
    ## convert sub-fields in 'x' to full R representation
{
    d <- suppressWarnings(as.integer(n))
    x <- if (!is.na(d)) {
        if (1L == d) x
        else {
            a <- array(unlist(strsplit(x, ",", fixed=TRUE)),
                       dim=c(d, nrow(x), ncol(x)),
                       dimnames=list(NULL, NULL, colnames(x)))
            aperm(a, c(2, 3, 1))
        }
    } else {
        split(strsplit(x, ",", fixed=TRUE), seq_len(nrow(x)))
    }
    switch(type,
           "Character"=,
           "String"=x,
           "Integer"={
               mode(x) <- "integer"
               x
           },
           "Float"={
               mode(x) <- "numeric"
               x
           },
           stop(sprintf("unhandled FORMAT type '%s'", type)))
}

.unpackVcfGeno <-
    function(geno, id, n, type)
    ## 'geno' is a named list,
    ## id, n, type are equal-length vectors of FORMAT info
{
    if (is.null(names(geno)))
        stop("'geno' must be a named list")
    Map(function(elt, nm, id, n, type) {
        if (is.na(idx <- match(nm, id))) {
            msg <- sprintf("element '%s' not found in FORMAT identifiers",
                           nm)
            stop(msg)
        }
        .unpackVcfField(elt, id[idx], n[idx], type[idx])
    }, geno, names(geno), MoreArgs=list(id, n, type))
}
