.STRAND_LEVELS <- c("+", "-", "*")

.onLoad <-
    function(libname, pkgname)
{
    if (!identical(levels(strand()), .STRAND_LEVELS))
        stop("internal: 'levels(strand())' not consistent with Rsamtools")
    .Call(.bamfile_init)
    .Call(.bcffile_init)
    .Call(.fafile_init)
    .Call(.tabixfile_init)
}
