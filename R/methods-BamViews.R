setMethod(BamViews, c(bamRanges="GRanges"), 
          function(bamPaths=character(0),
                   bamIndicies=bamPaths,
                   bamSamples=DataFrame(row.names=
                     make.unique(basename(bamPaths))),
                   bamRanges,
                   bamExperiment=list(), ...)
{
    new("BamViews", ..., bamPaths=bamPaths, bamIndicies=bamIndicies,
        bamSamples=bamSamples, bamRanges=bamRanges,
        bamExperiment=bamExperiment)
})

setMethod(BamViews, c(bamRanges="RangedData"),
          function(bamPaths=character(0),
                   bamIndicies=bamPaths,
                   bamSamples=DataFrame(row.names=
                     make.unique(basename(bamPaths))),
                   bamRanges,
                   bamExperiment=list(), ...)
{
    bamRanges <- as(bamRanges, "GRanges")
    callGeneric(bamPaths=bamPaths, bamIndicies=bamIndicies,
                bamSamples=bamSamples, bamRanges=bamRanges,
                bamExperiment=bamExperiment, ...)
})

setMethod(BamViews, c(bamRanges="missing"), 
          function(bamPaths=character(0),
                   bamIndicies=bamPaths,
                   bamSamples=DataFrame(row.names=
                     make.unique(basename(bamPaths))),
                   bamRanges,
                   bamExperiment=list(), ...,
                   auto.range=FALSE)
{
    if (length(bamPaths) != 0L && auto.range)
    {
        ## Guess ranges from BAM file headers
        pathsOk <- sapply(bamPaths, function(fl) {
            file.exists(fl) && !file.info(fl)$isdir
        })
        if (all(pathsOk)) {
            rngs <- lapply(scanBamHeader(bamPaths), "[[", "targets")
            nms <- unique(unlist(lapply(rngs, names), use.names=FALSE))
            ends <- sapply(nms, function(nm, rngs) {
                idx <- sapply(rngs, function(rng, nm) {
                    nm %in% names(rng)
                }, nm)
                if (sum(idx) > 0)
                    max(sapply(rngs[idx], "[[", nm))
                else
                    stop("Rsamtools internal: could not determine bamRanges")
            }, rngs)
            bamRanges <- GRanges(names(ends), IRanges(1L, ends))
        } else {
            warning("some files do not exist; bamRanges not defined")
            bamRanges <- GRanges()
        }
    } else {
        bamRanges <- GRanges()
    }
    callGeneric(bamPaths=bamPaths, bamIndicies=bamIndicies,
        bamSamples=bamSamples, bamRanges=bamRanges,
        bamExperiment=bamExperiment, ...)
})

setMethod(.validity, "BamViews", function(object) {
    msg <- NULL
    if (length(bamIndicies(object)) != length(bamPaths(object)))
        msg <- c(msg,
                 "length(bamIndicies(object)) != length(bamPaths(object))")
    if (length(bamPaths(object)) != nrow(bamSamples(object)))
        msg <- c(msg,
                 "length(bamPaths(object)) != nrow(bamSamples(object))")
    if (is.null(msg)) TRUE else msg
})

bamPaths <-
    function(x) slot(x, "bamPaths")

bamIndicies <-
    function(x) slot(x, "bamIndicies")

`bamDirname<-` <-
    function(x, ..., value)
{
    initialize(x,
               bamPaths=file.path(value, basename(bamPaths(x))),
               bamIndicies=file.path(value, basename(bamIndicies(x))))
}

bamSamples <-
    function(x) slot(x, "bamSamples")

`bamSamples<-` <-
    function(x, value) initialize(x, bamSamples=value)
        
bamRanges <-
    function(x) slot(x, "bamRanges")

`bamRanges<-` <-
    function(x, value) initialize(x, bamRanges=value)

bamExperiment <-
    function(x) slot(x, "bamExperiment")
    
setMethod(dim, "BamViews", function(x) {
    c(length(bamRanges(x)), length(bamPaths(x)))
})

setMethod(names, "BamViews", function(x) {
    rownames(bamSamples(x))
})

setReplaceMethod("names", "BamViews", function(x, value) {
    rownames(bamSamples(x)) <- value
    x
})

setMethod(dimnames, "BamViews", function(x) {
    list(names(bamRanges(x)),
         rownames(bamSamples(x)))
})

setReplaceMethod("dimnames", "BamViews", function(x, value) {
    names(bamRanges(x)) <- value[[1]]
    rownames(bamSamples(x)) <- value[[2]]
    x
})

setMethod("[", c("BamViews", "ANY", "missing"),
          function(x, i, j, ..., drop=TRUE)
{
    initialize(x, bamRanges=bamRanges(x)[i,])
})

setMethod("[", c("BamViews", "missing", "ANY"),
          function(x, i, j, ..., drop=TRUE)
{
    if (is.character(j))
        j <- match(j, colnames(x))
    if (any(is.na(j)))
        stop("subscript 'j' out of bounds")
    initialize(x, bamPaths=bamPaths(x)[j],
               bamIndicies=bamIndicies(x)[j],
               bamSamples=bamSamples(x)[j,,drop=FALSE])
})

setMethod("[", c("BamViews", "ANY", "ANY"),
          function(x, i, j, ..., drop=TRUE)
{
    if (is.character(i))
        j <- match(i, rownames(x))
    if (is.character(j))
        j <- match(j, colnames(x))
    if (any(is.na(i)))
        stop("subscript 'i' out of bounds")
    if (any(is.na(j)))
        stop("subscript 'j' out of bounds")
    initialize(x, bamRanges=bamRanges(x)[i,],
               bamPaths=bamPaths(x)[j],
               bamIndicies=bamIndicies(x)[j],
               bamSamples=bamSamples(x)[j,,drop=FALSE])
})

## action

## FIXME: temporary alias for srapply
.srapply <-
    function(X, FUN, ..., fapply=ShortRead:::.fapply(),
             reduce=ShortRead:::.reduce(), verbose=FALSE)
{
    require(ShortRead)
    result <- fapply(X, FUN, ..., verbose=verbose)
    reduce(result)
}

.BamViews_which <- function(file, param, missing)
{
    grange <- bamRanges(file)
    which <- split(ranges(grange), seqnames(grange))
    if (!missing && !identical(which, bamWhich(param)))
        warning("'bamRanges(file)' and 'bamWhich(param)' differ; using bamRanges(file)")
    which
}

.BamViews_delegate <-
    function(what, bamViews, fun, ...)
{
    idx <- structure(seq_len(ncol(bamViews)), names=names(bamViews))
    result <- .srapply(idx, fun, bamViews, ...)
    if (length(result) != ncol(bamViews)) {
        stop(sprintf("'%s' failed on '%s'", what,
                     paste(setdiff(names(bamViews), names(result)),
                                   collapse="' '")))
    }
    names(result) <- names(bamViews)
    do.call(new, list("SimpleList", listData=result,
                      elementMetadata=bamSamples(bamViews)))
}

setMethod(scanBam, "BamViews",
          function(file, index=file, ...,
                   param=ScanBamParam(what=scanBamWhat()))
{
    if (!missing(index))
        warning("using bamIndicies(file) for 'index'")
    bamWhich(param) <- .BamViews_which(file, param, missing(param))
    fun <- function(i, bamViews, ..., verbose)
        scanBam(file=bamPaths(bamViews)[i],
                index=bamIndicies(bamViews)[i], ...)
    .BamViews_delegate("scanBam", file, fun, ..., param=param)
})

setMethod(countBam, "BamViews",
          function(file, index=file, ..., param=ScanBamParam())
{
    if (!missing(index))
        warning("using bamIndicies(file) for 'index'")
    bamWhich(param) <- .BamViews_which(file, param, missing(param))
    fun <- function(i, bamViews, ..., verbose)
        countBam(file=bamPaths(bamViews)[i],
                 index=bamIndicies(bamViews)[i], ...)
    .BamViews_delegate("countBam", file, fun, ..., param=param)
})

setMethod(readBamGappedAlignments, "BamViews",
          function(file, index=file, use.names=FALSE, param=NULL)
{
    if (missing(index))
        index <- bamIndicies(file)
    if (is.null(param)) {
        param <- ScanBamParam(which=bamRanges(file))
    } else if (!identical(bamRanges(file), bamWhich(param))) {
        warning("'bamRanges(file)' and 'bamWhich(param)' differ; using 'bamRanges(file)'")
        bamWhich(param) <- bamRanges(file)
    }
    fun <- function(i, bamViews, verbose)
        readBamGappedAlignments(file=bamPaths(bamViews)[i],
                                index=bamIndicies(bamViews)[i],
                                use.names=use.names,
                                param=param)
    .BamViews_delegate("readBamGappedAlignments", file, fun)
})

## show
setMethod(show, "BamViews", function(object) {
    cat(class(object), "dim:",
        paste(dim(object), c("ranges", "samples"), collapse=" x "),
        "\n")
    cat("names:", IRanges:::selectSome(names(object)), "\n")
    cat("detail: use bamPaths(), bamSamples(), bamRanges(), ...",
        "\n")
})


## summarizeOverlaps methods

setMethod("summarizeOverlaps", c("BamViews", "missing"),
function(features, reads, mode, ignore.strand=FALSE, 
         ..., singleEnd=TRUE, param=ScanBamParam())
{
    se <- .processBamFiles(bamRanges(features), BamFileList(bamPaths(features)), 
        mode, ignore.strand, ..., singleEnd=singleEnd, param=param)
    nms <- rownames(colData(se))
    colData(se) <- DataFrame(colData(se), bamSamples(features),
        bamIndicies(features))
    rownames(colData(se)) <- nms
    exptData(se) <- SimpleList(bamExperiment=bamExperiment(features)) 
    se 
})
