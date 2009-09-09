setMethod(countBam, "character",
          function(file, index=file, ..., param=ScanBamParam())
{
    index <- normalizePath(path.expand(index))
    flag <- bamFlag(param)
    simpleCigar <- bamSimpleCigar(param)
    what <- bamWhat(param)
    which <- bamWhich(param)
    on.exit(.Call(.scan_bam_cleanup))
    if (!is.null(space(which))) {
        x <- .Call(.count_bam, file, index, "rbu",
                   list(space(which), start(which), end(which)),
                   flag, simpleCigar)
        data.frame(space=space(which), start=start(which),
                   end=end(which), width=width(which),
                   file=basename(file), records=x[["records"]],
                   nucleotides=x[["nucleotides"]])
     } else {
         x <- .Call(.count_bam, file, index, "rbu", NULL, flag,
                    simpleCigar)
         data.frame(space=NA, start=NA, end=NA, width=NA,
                   file=basename(file), records=x[["records"]],
                   nucleotides=x[["nucleotides"]])
     }
})
