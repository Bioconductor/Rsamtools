.STRAND_LEVELS <- c("+", "-", "*")

.onLoad <-
    function(libname, pkgname)
{
    if (!identical(levels(strand()), .STRAND_LEVELS))
        stop("internal: 'levels(strand())' not consistent with Rsamtools")

    msg <- "Rsamtools is experimental; expect frequent changes
            in data types and functionality.

            Issues with Windows remote-file retrieval exist; this
            is not yet reliable."

    message("\n",
            paste(strwrap(msg, indent=2, exdent=2), collapse="\n"),
            "\n")
}
