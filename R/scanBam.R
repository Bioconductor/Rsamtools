.scanBamTemplate <-
    function(tag=character(0))
{
    .Call(.scan_bam_template, tag)
}

.scanBam_template <-
    function(param)
{
    tmpl <- .scanBamTemplate(bamTag(param))
    tmpl[!names(tmpl) %in% c(bamWhat(param), "tag")] <- list(NULL)
    if (0L == length(tmpl[["tag"]]))
        tmpl["tag"] <- list(NULL)
    tmpl
}

.scanBam_postprocess <-
    function(x, param)
{
    which <- bamWhich(param)
    if (0L != length(space(which)))
        names(x) <-
            paste0(space(which), ":", .uunlist(start(which)), "-",
                   .uunlist(end(which)))
    lapply(x, Filter, f=Negate(is.null))
}

setMethod(scanBam, "character",
          function(file, index=file, ...,
                   param=ScanBamParam(what=scanBamWhat()))
{
    if (missing(index) && 0L == length(bamWhich(param)))
        index <- character(0)
    bam <- open(BamFile(file, index), "rb")
    on.exit(close(bam))
    callGeneric(bam, ..., param=param)
})
