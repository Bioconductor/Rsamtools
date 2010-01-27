BamViews <-
    function(bamPaths=character(0),
             bamSamples=new("DataFrame", nrows=length(bamPaths)),
             bamRanges=RangedData(),
             bamExperiment=list(), ...)
{
    new("BamViews",
        bamPaths=bamPaths,
        bamSamples=bamSamples,
        bamRanges=bamRanges,
        bamExperiment=bamExperiment)
}

setMethod(.validity, "BamViews", function(object) {
    msg <- NULL
    if (length(bamPaths(object)) != nrow(bamSamples(object)))
        msg <- c(msg,
                 "length(bamPaths(object)) != nrow(bamSamples(object))")
    if (is.null(msg)) TRUE else msg
})

bamPaths <-
    function(object) slot(object, "bamPaths")

bamSamples <-
    function(object) slot(object, "bamSamples")

bamRanges <-
    function(object) slot(object, "bamRanges")

bamExperiment <-
    function(object) slot(object, "bamExperiment")
    
setMethod(dim, "BamViews", function(x) {
    c(nrow(bamRanges(x)), length(bamPaths(x)))
})

setMethod("[", c("BamViews", "ANY", "missing"),
          function(x, i, j, ..., drop=TRUE)
{
    initialize(x, bamRanges=bamRanges(x)[i,])
})

setMethod("[", c("BamViews", "missing", "ANY"),
          function(x, i, j, ..., drop=TRUE)
{
    initialize(x, bamPaths=bamPaths(x)[j],
               bamSamples=bamSamples(x)[j,,drop=FALSE])
})

setMethod("[", c("BamViews", "ANY", "ANY"),
          function(x, i, j, ..., drop=TRUE)
{
    initialize(x, bamRanges=bamRanges(x)[i,],
               bamPaths=bamPaths(x)[j],
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
        olaps <- findOverlaps(bamRanges(x), rangesList, ...)
        i <- queryHits(olaps)
        if (length(i)==0L)
            return(RangedData())
        rng <- ranges(olaps, ranges(bamRanges(x)), rangesList)
        if (sum(sapply(rng, length)) != length(i))
            stop("Rsamtools internal: 'ranges()' and 'findOverlaps' return different numbers of matches")
        RangedData(rng, values(bamRanges(x)[i,]))
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
    initialize(x,
               bamRanges=.subset_RangedData(x, i, ...),
               bamPaths=bamPaths(x)[j],
               bamSamples=bamSamples(x)[j,])
})

setMethod(show, "BamViews", function(object) {
    cat(class(object), "dim:",
        paste(dim(object), c("ranges", "samples"), collapse=" x "),
        "\n")
    cat("bamPaths:",
        if (length(bamPaths(object)) != 0) "\n  ",
        paste(IRanges:::selectSome(bamPaths(object)), collapse="\n  "),
        "\n", sep="")
    cat("detail: use bamSamples(), bamRanges(), bamExperiment()",
        "\n")
})
