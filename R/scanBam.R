.scanBamTemplate <- function()
    .Call(.scan_bam_template)

setMethod(scanBam, "character",
          function(file, index=file, ..., param=ScanBamParam())
{
    tmpl <- .scanBamTemplate()
    tmpl[!names(tmpl) %in% bamWhat(param)] <- list(NULL)
    x <- .io_bam(.scan_bam, file, index, tmpl, param=param)
    which <- bamWhich(param)
    if (!is.null(space(which)))
        names(x) <-
            paste(space(which), ":", start(which), "-", end(which),
                  sep="")
    x
})
