.scanBamTemplate <-
    function(seqlevels=factor(), tag=character(0))
{
    .Call(.scan_bam_template, seqlevels, tag)
}

.scanBam_template <-
    function(file, param)
{
    tmpl <- .scanBamTemplate(factor(levels=seqlevels(file)), bamTag(param))
    ## set those elements of the template that are not 'tag' (treat
    ## specially because nested list) nor specified by 'what'
    ## parameter to NULL
    tmpl[!names(tmpl) %in% c(bamWhat(param), "tag")] <- list(NULL)
    if (0L == length(tmpl[["tag"]]))
        tmpl["tag"] <- list(NULL)
    tmpl
}

## return rname:start-end values for outer list elements
.scanBam_extract_which_labels <- function(param) {
    which <- bamWhich(param)
    if( 0L != length(space(which)))
        paste0(space(which), ":", .uunlist(start(which)), "-",
               .uunlist(end(which)))
    else
        NULL
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
    scanBam(bam, ..., param=param)
})
