setMethod("testPairedEndBam", "character",
    function(file, index=file, ...)
{
    bf <- BamFile(file, index)
    testPairedEndBam(bf, ...)
})

setMethod("testPairedEndBam", "BamFile",
    function(file, index=file, ...)
{
    yieldSize <- yieldSize(file)
    yieldSize(file) <- 1000000L
    on.exit(yieldSize(file) <- yieldSize)
    isPaired <- FALSE
    repeat {
        flag <- scanBam(file, param=ScanBamParam(what="flag"))[[1]]$flag
        isPaired <- any(bamFlagTest(flag, "isPaired"))
        if (isPaired || length(flag) == 0)
            break
    }
    isPaired
})
