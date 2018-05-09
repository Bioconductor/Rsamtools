.io_bam <-
    function(func, file, ..., param)
{
    flag <- bamFlag(param, asInteger=TRUE)
    simpleCigar <- bamSimpleCigar(param)
    tagFilter <- bamTagFilter(param)
    mapqFilter <- bamMapqFilter(param)
    which <- bamWhich(param)
    if (!all(names(which) %in% seqlevels(file))) {
        bad <- setdiff(names(which), seqlevels(file))
        stop("seqlevels(param) not in BAM header:",
             "\n    seqlevels: ", paste(sQuote(bad), collapse=", "),
             "\n    file: ", path(file),
             "\n    index: ", index(file))
    }
    regions <-
        if (0L != length(space(which)))
            list(as.character(space(which)), .uunlist(start(which)),
                 .uunlist(end(which)))
        else NULL
    on.exit(.Call(.scan_bam_cleanup))

    .io_check_exists(path(file))
    tryCatch({
        .Call(func, .extptr(file), regions, flag, simpleCigar, tagFilter,
              mapqFilter, ...)
    }, error=function(err) {
        stop(conditionMessage(err), "\n  file: ", path(file),
             "\n  index: ", index(file))
    })
}
