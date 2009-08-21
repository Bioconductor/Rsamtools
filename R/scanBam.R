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
    ## file
    bam <- .Call(.scan_bam_open, file, "rbu")
    ## which
    if (!is.null(space(which))) {
        result <- mapply(function(bam, tmpl, space, start, end) {
            spc <- list(space, start, end)
            on.exit(.Call(.scan_bam_cleanup))
            .Call(.scan_bam, bam, tmpl, spc, flag, simpleCigar)
        }, space(which), start(which), end(which),
                         MoreArgs=list(bam=bam, tmpl=tmpl),
                         USE.NAMES=FALSE,
                         SIMPLIFY=FALSE)
        names(result) <-
            paste(space(which), ":", start(which), "-", end(which),
                  sep="")
        result
     } else {
         on.exit(.Call(.scan_bam_cleanup))
         list(.Call(.scan_bam, bam, tmpl, NULL, flag, simpleCigar))
     }
})
