setMethod(scanBcfHeader, "character",
    function(file, ...)
{
    Map(function(file, mode) {
        bf <- open(BcfFile(file, character(0), ...))
        on.exit(close(bf))
        scanBcfHeader(bf)
    }, file, ...)
})
