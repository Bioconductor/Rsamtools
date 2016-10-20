.normalizeRangesList <-
    function(rangesList)
{
    nms <- names(rangesList)
    reducedList <- if (0 != length(rangesList) && is.null(nms)) {
        ## special case, all names missing
        rng <- Reduce(append, as(rangesList, "list"))
        IRangesList(reduce(rng, drop.empty.ranges=TRUE))
    } else if (any(duplicated(nms))) {
        unms <- unique(nms)
        lst <- lapply(unms, function(nm, rnglist) {
            idx <- names(rnglist) == nm
            rng <- Reduce(append, as(rnglist[idx], "list"))
            reduce(rng, drop.empty.ranges=TRUE)
        }, rnglist=rangesList)
        names(lst) <- unms
        do.call(IRangesList, lst)
    } else {
        reduce(rangesList, drop.empty.ranges=TRUE)
    }
    reducedList[lengths(reducedList) != 0]
}

.filterBam_preprocess <-
    function(file, param)
{
    which <- .normalizeRangesList(bamWhich(param))
    hnames <- seqlevels(file)
    o <- order(match(names(which), hnames))
    what <- bamWhat(param)
    if (asMates(file))
        what <- union(what, c("mate_status", "groupid"))
    initialize(param, which=which[o], what=what)
}

setMethod(filterBam, "character",
          function(file, destination, index=file, ...,
                   filter=FilterRules(),
                   indexDestination=TRUE,
                   param=ScanBamParam(what=scanBamWhat()))
{
    if (missing(index) && 0L == length(bamWhich(param)))
        index <- character(0)
    bam <- open(BamFile(file, index), "rb")
    on.exit(close(bam))
    filterBam(bam, destination, ..., filter=filter,
              indexDestination=indexDestination, param=param)
})
