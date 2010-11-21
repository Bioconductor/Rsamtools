.normalizeRangesList <-
    function(rangesList)
{
    nms <- names(rangesList)
    if (0 != length(rangesList) && is.null(nms)) {
        ## special case, all names missing
        rng <- as(Reduce(append, as(rangesList, "list")),
                  "NormalIRanges")
        RangesList(rng)
    } else if (any(duplicated(nms))) {
        unms <- unique(nms)
        lst <- lapply(unms, function(nm, rnglist) {
            idx <- names(rnglist) == nm
            as(Reduce(append, as(rnglist[idx], "list")),
               "NormalIRanges")
        }, rnglist=rangesList)
        names(lst) <- unms
        do.call(RangesList, lst)
    } else {
        endoapply(rangesList, as, "NormalIRanges")
    }
}

.filterBam_preprocess <-
    function(file, param)
{
    which <- .normalizeRangesList(bamWhich(param))
    hnames <- names(scanBamHeader(file)[[1]][["length"]])
    o <- order(match(names(which), hnames))
    initialize(param, which=which[o])
}

setMethod(filterBam, "character",
          function(file, destination, index=file, ...,
                   indexDestination=TRUE, param=ScanBamParam())
{
    index <- 
        if (missing(index) && 0L == length(bamWhich(param)))
            character(0)
        else .normalizePath(index)
    bam <- openBam(.normalizePath(file), index, "rb")
    on.exit(close(bam))
    callGeneric(bam, destination, ...,
                indexDestination=indexDestination, param=param)
})
