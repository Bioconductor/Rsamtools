BamViews <-
    function(bamPaths=character(0),
             bamIndicies=bamPaths,
             bamSamples=new("DataFrame", nrows=length(bamPaths),
                            rownames=make.unique(basename(bamPaths))),
             bamRanges=RangedData(),
             bamExperiment=list(), ...,
             auto.range=FALSE)
{
    if (length(bamPaths) != 0L && missing(bamRanges) && auto.range)
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
            bamRanges <-
                RangedData(IRanges(1L, ends), space=names(ends))
        } else {
            warning("some files do not exist; bamRanges not defined")
        }
    }
    new("BamViews", ..., bamPaths=bamPaths, bamIndicies=bamIndicies,
        bamSamples=bamSamples, bamRanges=bamRanges,
        bamExperiment=bamExperiment)
}

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
    c(nrow(bamRanges(x)), length(bamPaths(x)))
})

setMethod(names, "BamViews", function(x) {
    rownames(bamSamples(x))
})

setReplaceMethod("names", "BamViews", function(x, value) {
    rownames(bamSamples(x)) <- value
    x
})

setMethod(dimnames, "BamViews", function(x) {
    list(rownames(bamRanges(x)),
         rownames(bamSamples(x)))
})

setReplaceMethod("dimnames", "BamViews", function(x, value) {
    rownames(bamRanges(x)) <- value[[1]]
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

.subset_RangedData <-
    function(x, rangesList, ...)
{
    stop("'[,BamViews,RangesList-method' not yet supported")
    ## FIXME: nrow(RangesList())==0 is a proxy for 'all ranges'
    if (length(rangesList) == 0L)
        bamRanges(x)
    else if (nrow(bamRanges(x)) == 0L)
        RangedData(rangesList)
    else {
        bamRanges(x)[bamRanges(x) %in% rangesList,]
    }
}

setMethod("[", c("BamViews", "RangesList", "missing"),
          function(x, i, j, ..., drop=TRUE)
{
    initialize(x, bamRanges=.subset_RangedData(x, i, ...))
})

setMethod("[", c("BamViews", "RangesList", "ANY"),
          function(x, i, j, ..., drop=TRUE)
{
    if (is.character(j))
        j <- match(j, colnames(x))
    if (any(is.na(j)))
        stop("subscript 'j' out of bounds")
    initialize(x,
               bamRanges=.subset_RangedData(x, i, ...),
               bamPaths=bamPaths(x)[j],
               bamSamples=bamSamples(x)[j,])
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

setMethod(readBAMasGappedAlignments, "BamViews",
          function(file, index, ..., which=RangesList(),
                   ans.subtype="Alignments0")
{
    if (missing(index))
        index <- bamIndicies(file)
    if (!missing(which))
        file <- file[which,]

    fun <- function(i, bamViews, ..., verbose)
        readBAMasGappedAlignments(file=bamPaths(bamViews)[i],
                                  index=bamIndicies(bamViews)[i],
                                  ...,
                                  ans.subtype=ans.subtype)
    idx <- structure(seq_len(ncol(file)), names=names(file))
    res <- .srapply(idx, fun, bamViews=file, ...,
                    which=ranges(bamRanges(file)))
    if (length(res) != ncol(file))
        stop("'readBAMasGappedAlignments' failed on '",
             paste(setdiff(names(file), names(res)), collapse="' '"),
             "'")

    names(res) <- rownames(bamSamples(file))
    do.call(new, list("SimpleList", listData=res,
                      elementMetadata=bamSamples(file)))
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
