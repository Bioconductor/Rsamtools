.BamSampler <- setRefClass("BamSampler", contains="BamFile")

BamSampler <-
    function (file, index = file, ..., yieldSize, obeyQname = FALSE,
              asMates = FALSE) 
{
    .RsamtoolsFile(.BamSampler, .normalizePath(file), .normalizePath(index), 
        yieldSize = yieldSize, obeyQname = obeyQname, asMates = asMates, ...)
}

setMethod("scanBam", "BamSampler",
    function(file, index=file, ...,
             param=ScanBamParam(what=scanBamWhat()))
{
    tot <- sampleSize <- yieldSize(file)
    if (is.na(yieldSize(file)))
        stop("'yieldSize' must not be NA")

    bfile <- as(file, "BamFile")
    open(bfile, "rb")
    on.exit(close(bfile))

    smpl <- .quickUnlist(scanBam(bfile, param=param))
    repeat {
        yld <- .quickUnlist(scanBam(bfile, param=param))
        yld_n <- length(yld[[1]])
        if (length(yld[[1]]) == 0L)
            break
        tot <- tot + yld_n
        keep <- rbinom(1L, yld_n, yld_n/ tot)
        if (keep == 0L)
            next

        i <- sample(sampleSize, keep)
        j <- sample(yld_n, keep)
        smpl <- Map(function(x, y, i, j) {
            x[i] <- y[j]
            x
        }, smpl, yld, MoreArgs=list(i=i, j=j))
    }
    list(smpl)
})

setMethod(show, "BamSampler", function(object) {
    callNextMethod()
})
