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
    if (is.na(yieldSize))
        yieldSize(file) <- 1000000L
    on.exit(yieldSize(file) <- yieldSize)
    if (isOpen(file))
        close(file)
    open(file)
    isPaired <- FALSE
    tot <- 0
    repeat {
        flag <- scanBam(file, param=ScanBamParam(what="flag"))[[1]]$flag
        isPaired <- any(bamFlagTest(flag, "isPaired"))
        if (isPaired || length(flag) == 0L)
            break
        tot <- tot + length(flag)
        message(tot, " ", appendLF=FALSE)
    }
    isPaired
})
