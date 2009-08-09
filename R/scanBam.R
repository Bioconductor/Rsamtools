.scanBamTemplate <- function()
    .Call(ShortRead:::.scan_bam_template)

setMethod(scanBam, "character",
          function(file, ..., param=ScanBamParam())
{
    flag <- bamFlag(param)
    what <- bamWhat(param)
    which <- bamWhich(param)
    ## template
    tmpl <- .scanBamTemplate()
    if (!all(what %in% names(tmpl)))
        warn("'what' argument contains invalid names:\n  ",
             paste(what[!what %in% names(tmpl)], collapse=", "))
    tmpl[!names(tmpl) %in% what] <- list(NULL)
    ## file
    bam <- .Call(ShortRead:::.scan_bam_open, file, "rbu")
    ## which
    if (!is.null(space(which))) {
        result <- mapply(function(bam, tmpl, space, start, end) {
            spc <- list(space, start, end)
            .Call(ShortRead:::.scan_bam, bam, tmpl, spc, flag)
        }, space(which), start(which), end(which),
                         MoreArgs=list(bam=bam, tmpl=tmpl),
                         USE.NAMES=FALSE,
                         SIMPLIFY=FALSE)
        names(result) <-
            paste(space(which), ":", start(which), "-", end(which),
                  sep="")
        result
     } else 
        .Call(ShortRead:::.scan_bam, bam, tmpl, NULL, flag)
})
