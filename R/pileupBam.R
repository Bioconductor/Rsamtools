.PileupParam <-
    function(flag=scanBamFlag(),
             minBaseQuality=13L,
             minMapQuality=0L,
             minDepth=0L,
             maxDepth=250L,
             yieldSize=1L,
             yieldBy=c("range", "position"),
             yieldAll=FALSE,
             which=GRanges(),
             what=c("seq", "qual"))
{
    yieldBy <- match.arg(yieldBy)
    if ("range" == yieldBy && yieldSize != 1)
        stop("'yieldSize' must equal 1 when 'yieldBy=\"range\"'")
    as.list(environment())
}

setGeneric(".asSpace", function(x) standardGeneric(".asSpace"))

setMethod(.asSpace, "RangesList", function(x) {
    list(as.character(space(x)), .uunlist(start(x)), .uunlist(end(x)))
})

setMethod(.asSpace, "GRanges", function(x) {
    list(as.character(seqnames(x)), start(x), end(x))
})

.pileupBam <-
    function(files, FUN=identity, ..., param)
{
    FUN <- match.fun(FUN)
    tryCatch({
        files <- lapply(files, .extptr)
        param <- as.list(param)
        space <-
            if (0L != length(param[["which"]])) .asSpace(param[["which"]])
            else NULL
        what <- logical(2);
        param[["what"]] <- c("seq", "qual") %in% param[["what"]]
        .Call(.pileup_bam, files, space, param, FUN)
    }, error=function(err) {
        stop("pileupBam: ", conditionMessage(err), call.=FALSE)
    })
}

## setGeneric("pileupBam",
##     function(file, index=file, ..., param=PileupParam())
##         standardGeneric("pileupBam"),
##     signature="file")
