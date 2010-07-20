.io_bam <-
    function(func, file, index, ..., param)
{
    idx <- !grepl("^(ftp)|(http)://", file)
    if (!all(sapply(file[idx], file.exists))) {
        msg <- paste(sprintf("'%s'", file[idx]), collapse="\n  ")
        stop("'file' elements do not exist:\n  ", msg)
    }
    file <- .normalizePath(file)
    index <- .normalizePath(index)
    flag <- bamFlag(param)
    simpleCigar <- bamSimpleCigar(param)
    what <- bamWhat(param)
    tmpl <- .scanBamTemplate()
    if (!all(what %in% names(tmpl)))
        warning("'what' argument contains invalid names:\n  ",
             paste(what[!what %in% names(tmpl)], collapse=", "))
    which <- bamWhich(param)
    on.exit(.Call(.scan_bam_cleanup))
    if (0L != length((space(which))))
        .Call(func, file, index, "rb",
              list(as.character(space(which)), .uunlist(start(which)),
                   .uunlist(end(which))),
              flag, simpleCigar, ...)
    else 
        .Call(func, file, index, "rb", NULL, flag,
              simpleCigar, ...)
}
