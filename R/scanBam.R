.scanBamTemplate <-
    function(tag=character(0))
{
    .Call(.scan_bam_template, tag)
}

setMethod(scanBam, "character",
          function(file, index=file, ..., param=ScanBamParam())
{
    tmpl <- .scanBamTemplate(bamTag(param))
    tmpl[!names(tmpl) %in% c(bamWhat(param), "tag")] <- list(NULL)
    if (0L == length(tmpl[["tag"]]))
        tmpl["tag"] <- list(NULL)
    reverseComplement <- bamReverseComplement(param)
    x <- .io_bam(.scan_bam, file, index, reverseComplement, tmpl,
                 param=param)
    which <- bamWhich(param)
    if (!is.null(space(which)))
        names(x) <-
            paste(space(which), ":", .uunlist(start(which)), "-",
                  .uunlist(end(which)), sep="")
    lapply(x, Filter, f=Negate(is.null))
})
