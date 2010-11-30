setMethod(scanBcf, "character",
    function(file, index=file, ..., param=ScanBcfParam())
{
    bcf <- open(BcfFile(file, index))
    on.exit(close(bcf))
    callGeneric(bcf, ..., param=param)
})
