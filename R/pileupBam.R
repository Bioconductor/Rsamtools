.PileupParam <-
    function(flag=scanBamFlag(),
             minBaseQuality=13L,
             minMapQuality=0L,
             minDepth=0L,
             maxDepth=250L,
             which=GRanges())
{
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
    function(files, callback=identity, ..., param)
{
    tryCatch({
        files <- lapply(files, .extptr)
        param <- as.list(param)
        space <-
            if (0L != length(param[["which"]])) .asSpace(param[["which"]])
            else NULL
        .Call(.pileup_bam, files, space, param, callback)
    }, error=function(err) {
        stop("pileupBam: ", conditionMessage(err), call.=FALSE)
    })
}

## setGeneric("pileupBam",
##     function(file, index=file, ..., param=PileupParam())
##         standardGeneric("pileupBam"),
##     signature="file")
