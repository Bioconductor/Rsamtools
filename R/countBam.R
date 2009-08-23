setMethod(countBam, "character",
          function(file, ..., param=ScanBamParam())
{
    flag <- bamFlag(param)
    simpleCigar <- bamSimpleCigar(param)
    what <- bamWhat(param)
    which <- bamWhich(param)
    bam <- .Call(.scan_bam_open, file, "rbu")
    ## FIXME: error handling
    if (!is.null(space(which))) {
        result <- mapply(function(bam, space, start, end) {
            spc <- list(space, start, end)
            .Call(.count_bam, bam, spc, flag, simpleCigar)
        }, space(which), start(which), end(which),
                         MoreArgs=list(bam=bam),
                         USE.NAMES=FALSE)
        names(result) <-
            paste(space(which), ":", start(which), "-", end(which),
                  sep="")
        result
     } else {
         .Call(.count_bam, bam, NULL, flag, simpleCigar)
     }
})
