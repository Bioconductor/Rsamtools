## internal function; not robust enough for public use
##
## determine format (SAM / BAM / CRAM / VCF / BCF ...) of file
##
## @param x character(1) file path.
.fileFormat <-
    function(x)
{
    stopifnot(is.character(x), length(x) == 1L, nzchar(x))
    .Call(.file_format, x)
}
