.io_bam <-
    function(func, file, ..., param)
{
    flag <- bamFlag(param, asInteger=TRUE)
    simpleCigar <- bamSimpleCigar(param)
    tagFilter <- bamTagFilter(param)
    which <- bamWhich(param)
    space <- 
        if (0L != length(space(which)))
            list(as.character(space(which)), .uunlist(start(which)),
                 .uunlist(end(which)))
        else NULL
    on.exit(.Call(.scan_bam_cleanup))

    .io_check_exists(path(file))
    tryCatch({
        .Call(func, .extptr(file), space, flag, simpleCigar, tagFilter, ...)
    }, error=function(err) {
        stop(conditionMessage(err), "\n  file: ", path(file),
             "\n  index: ", index(file))
    })
}
