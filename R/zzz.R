.onLoad <-
    function(libname, pkgname)
{
    msg <- "Rsamtools is experimental; expect frequent changes
            in data types and functionality"
    message("\n",
            paste(strwrap(msg, indent=2, exdent=2), collapse="\n"),
            "\n")
}
