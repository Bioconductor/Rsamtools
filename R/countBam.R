.countBam_postprocess <- function(x, file, param)
{
    which <- bamWhich(param)
    bfile <- basename(path(file))
    if (0L != length(space(which))) {
        data.frame(space=space(which), start=.uunlist(start(which)),
                   end=.uunlist(end(which)),
                   width=.uunlist(width(which)), file=bfile,
                   records=x[["records"]],
                   nucleotides=x[["nucleotides"]])
    } else {
        data.frame(space=NA, start=NA, end=NA, width=NA,
                   file=bfile, records=x[["records"]],
                   nucleotides=x[["nucleotides"]])
    }
}

setMethod(countBam, "character",
          function(file, index=file, ..., param=ScanBamParam())
{
    index <- 
        if (missing(index) && 0L == length(bamWhich(param)))
            character(0)
        else .normalizePath(index)
    bam <- open(BamFile(file, index), "rb")
    on.exit(close(bam))
    callGeneric(bam, ..., param=param)
})
