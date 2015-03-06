.BamSampler <- setRefClass("BamSampler", contains="BamFile")

BamSampler <-
    function (file, index = file, ..., yieldSize, obeyQname = FALSE,
              asMates = FALSE, qnamePrefixEnd = NA,
              qnameSuffixStart = NA) 
{
    .Deprecated(msg=paste0("'BamSampler' is deprecated. Use 'REDUCEsampler' ",
                           "in the GenomicFiles package"))
    qnamePrefixEnd <- .check_qname_arg(qnamePrefixEnd, "qnamePrefixEnd")
    qnameSuffixStart <- .check_qname_arg(qnameSuffixStart, "qnameSuffixStart")
    .RsamtoolsFile(.BamSampler, .normalizePath(file), .normalizePath(index), 
        yieldSize = yieldSize, obeyQname = obeyQname, asMates = asMates, 
        qnamePrefixEnd = qnamePrefixEnd, qnameSuffixStart = qnameSuffixStart, 
        ...)
}

setMethod("scanBam", "BamSampler",
    function(file, index=file, ...,
             param=ScanBamParam(what=scanBamWhat()))
{
    if (0L == length(bamWhat(param)) && 0L == length(bamTag(param))) {
        txt <- "no BAM fields selected for input (niether 'bamWhat(param)'
                nor 'bamTag(param)' defined)"
        stop(paste(strwrap(txt), collapse="\n  "))
    }

    sampleSize <- yieldSize(file)
    if (is.na(yieldSize(file)))
        stop("'yieldSize' must not be NA")

    bfile <- as(file, "BamFile")
    open(bfile, "rb")
    on.exit(close(bfile))

    smpl <- S4Vectors:::quick_unlist(unname(scanBam(bfile, param=param)))
    tot <- length(smpl[[1]])
    if (tot > sampleSize) {             # e.g., ranges
        idx <- sample(tot, sampleSize)
        smpl <- lapply(smpl, `[`, idx)
    }
    repeat {
        yld <- S4Vectors:::quick_unlist(scanBam(bfile, param=param))
        yld_n <- length(yld[[1]])
        if (length(yld[[1]]) == 0L)
            break
        tot <- tot + yld_n
        keep <- rbinom(1L, yld_n, yld_n / tot)
        if (keep == 0L)
            next

        i <- sample(sampleSize, keep)
        j <- sample(yld_n, keep)
        smpl <- Map(function(x, y, i, j) {
            x[i] <- y[j]
            x
        }, smpl, yld, MoreArgs=list(i=i, j=j))
    }
    lst <- list(smpl)
    attr(lst, "BamSamplerStatistics") <-
        c(yieldSize=sampleSize, totalRead=tot, yield=length(smpl[[1]]))
    lst
})

setMethod(show, "BamSampler", function(object) {
    callNextMethod()
})
