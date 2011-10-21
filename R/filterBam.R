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
        if (is(rangesList, "CompressedIRangesList"))
            rangesList <- as(rangesList, "SimpleIRangesList")
        endoapply(rangesList, as, "NormalIRanges")
    }
}

.filterBam_preprocess <-
    function(file, param)
{
    which <- .normalizeRangesList(bamWhich(param))
    hnames <- names(scanBamHeader(file)[["targets"]])
    o <- order(match(names(which), hnames))
    initialize(param, which=which[o])
}

setMethod(filterBam, "character",
          function(file, destination, index=file, ...,
                   indexDestination=TRUE,
                   param=ScanBamParam(what=scanBamWhat()))
{
    if (missing(index) && 0L == length(bamWhich(param)))
        index <- character(0)
    bam <- open(BamFile(file, index), "rb")
    on.exit(close(bam))
    callGeneric(bam, destination, ...,
                indexDestination=indexDestination, param=param)
})
