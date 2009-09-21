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

setMethod(filterBam, "character",
          function(file, destination, index=file, ...,
                   indexDestination=TRUE, param=ScanBamParam())
{
    if (1 != length(file) || !is.character(file))
        stop("'file' must be character(1)")
    which <- .normalizeRangesList(bamWhich(param))
    hnames <- names(readBamHeader(file)[[1]][["length"]])
    o <- order(match(names(which), hnames))
    param <- initialize(param, which=which[o])
    destination <- .normalizePath(destination)
    mode <- "wb"
    fl <- .io_bam(.filter_bam, file, index, destination, mode,
                  param=param)
    if (indexDestination)
        indexBam(fl)
    fl
})
