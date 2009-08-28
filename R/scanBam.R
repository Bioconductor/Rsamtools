.scanBamTemplate <- function()
    .Call(.scan_bam_template)

setMethod(scanBam, "character",
          function(file, ..., param=ScanBamParam())
{
    flag <- bamFlag(param)
    simpleCigar <- bamSimpleCigar(param)
    what <- bamWhat(param)
    which <- bamWhich(param)
    ## template
    tmpl <- .scanBamTemplate()
    if (!all(what %in% names(tmpl)))
        warn("'what' argument contains invalid names:\n  ",
             paste(what[!what %in% names(tmpl)], collapse=", "))
    tmpl[!names(tmpl) %in% what] <- list(NULL)
    ## which
    on.exit(.Call(.scan_bam_cleanup))
    if (!is.null(space(which))) {
        result <- .Call(.scan_bam, file, "rbu", tmpl,
                        list(space(which), start(which), end(which)),
                        flag, simpleCigar)
        names(result) <-
            paste(space(which), ":", start(which), "-", end(which),
                  sep="")
        result
     } else {
         .Call(.scan_bam, file, "rbu", tmpl, NULL, flag, simpleCigar)
     }
})
