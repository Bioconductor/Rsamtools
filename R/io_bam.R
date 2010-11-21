.io_bam_check_exists <-
    function(file)
{
    idx <- !grepl("^(ftp)|(http)://", file)
    if (!all(sapply(file[idx], file.exists))) {
        msg <- paste(sprintf("'%s'", file[idx]), collapse="\n  ")
        stop("'file' elements do not exist:\n  ", msg)
    }
}

.io_bam <-
    function(func, file, ..., param)
{
    .io_bam_check_exists(bamPath(file))
    flag <- bamFlag(param)
    simpleCigar <- bamSimpleCigar(param)
    which <- bamWhich(param)
    space <- 
        if (0L != length((space(which))))
            list(as.character(space(which)), .uunlist(start(which)),
                 .uunlist(end(which)))
        else NULL
    on.exit(.Call(.scan_bam_cleanup))
    tryCatch({
        .Call(func, .extptr(file), space, flag, simpleCigar, ...)
    }, error=function(err) {
        stop(conditionMessage(err), "\n  file: ", bamPath(file))
    })
}
