setMethod("testPairedEndBam", "character",
    function(file, index=file, ..., fast=TRUE)
{
    bf <- BamFile(file, index)
    testPairedEndBam(bf, ..., fast=fast)
})

setMethod("testPairedEndBam", "BamFile",
    function(file, index=file, ..., fast=TRUE)
{
    yieldSize(file) <- if (fast) 1000000 else NA_integer_
    flag <- scanBam(file, param=ScanBamParam(what="flag"))[[1]]$flag
    any(bamFlagTest(flag, "isPaired"))
})