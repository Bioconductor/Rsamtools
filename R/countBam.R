setMethod(countBam, "character",
          function(file, index=file, ..., param=ScanBamParam())
{
    x <- .io_bam(.count_bam, file, index, param=param)
    which <- bamWhich(param)
    if (!is.null(space(which))) {
        data.frame(space=space(which), start=start(which),
                   end=end(which), width=width(which),
                   file=basename(file), records=x[["records"]],
                   nucleotides=x[["nucleotides"]])
    } else {
        data.frame(space=NA, start=NA, end=NA, width=NA,
                   file=basename(file), records=x[["records"]],
                   nucleotides=x[["nucleotides"]])
    }
})
