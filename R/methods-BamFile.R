setMethod(isOpen, "BamFile",
    function(con, rw="")
{
    if (!missing(rw) && rw == "read")
        stop("'rw' must be 'read'")
    .Call(.bamfile_isopen, .extptr(con))
})

BamFile <-
    function(file, index=file)
{
    .RsamtoolsFile(.BamFile, file, index)
}

open.BamFile <-
    function(con, ...)
{
    .io_bam_check_exists(path(con))
    con$.extptr <- .Call(.bamfile_open, path(con), index(con), "rb")
    invisible(con)
}

close.BamFile <-
    function(con, ...)
{
    if (!isOpen(con))
        stop("isOpen(<BamFile>) is not 'TRUE'")
    con$.extptr <- .Call(.bamfile_close, .extptr(con))
    invisible(con)
}

## scanBam, filterBam, countBam

setMethod(scanBamHeader, "BamFile",
          function(files, ...)
{
    .Call(.read_bamfile_header, .extptr(files))
})

setMethod(scanBam, "BamFile",
          function(file, index=file, ..., param=ScanBamParam())
{
    if (!isOpen(file))
        stop("BamFile not open")
    if (!missing(index))
        warning("'index' ignored for scanBam,BamFile-method")
    reverseComplement <- bamReverseComplement(param)
    tmpl <- .scanBam_template(param)
    x <- .io_bam(.scan_bamfile, file, param=param,
                 path(file), index(file), "rb",
                 reverseComplement, tmpl)
    .scanBam_postprocess(x, param)
})

setMethod(countBam, "BamFile",
          function (file, index=file, ..., param = ScanBamParam())
{
    if (!isOpen(file))
        stop("BamFile not open")
    if (!missing(index))
        warning("'index' ignored for countBam,BamFile-method")
    x <- .io_bam(.count_bamfile, file, param=param)
    .countBam_postprocess(x, file, param)
})

setMethod(filterBam, "BamFile",
          function (file, destination, index=file, ...,
                    indexDestination=TRUE, param=ScanBamParam())
{
    if (!isOpen(file))
        stop("BamFile not open")
    if (!missing(index))
        warning("'index' ignored for filterBam,BamFile-method")
    param <- .filterBam_preprocess(file, param)
    destination <- .normalizePath(destination)
    fl <- .io_bam(.filter_bamfile, file, param=param, destination, "wb")
    if (indexDestination)
        indexBam(fl)
    fl
})

setMethod(indexBam, "BamFile", function(files, ...) {
    callGeneric(path(files), ...)
})

setMethod(sortBam, "BamFile",
    function(file, destination, ..., byQname=FALSE, maxMemory=512)
{
    callGeneric(path(files), destination, ...,
                byQname=byQname, maxMemory=maxMemory)
})


setMethod(readBamGappedAlignments, "BamFile",
          function(file, index, ..., which=RangesList())
{
    param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE,
                            isDuplicate=FALSE),
                          what=c("rname", "strand", "pos", "cigar"),
                          which=which)
    bam <- scanBam(file, param=param)
    ## unlist(list(factor())) returns integer(0), so exit early if all
    ## values are empty
    bam <- bam[sapply(bam, function(x) length(x$rname) != 0)]
    if (0L == length(bam))
        return(GappedAlignments())
    rname <- unlist(unname(lapply(bam, "[[", "rname")))
    strand <- unlist(unname(lapply(bam, "[[", "strand")))
    pos <- unlist(unname(lapply(bam, "[[", "pos")))
    cigar <- unlist(unname(lapply(bam, "[[", "cigar")))
    ## Calls the appropriate constructor.
    GappedAlignments(rname=rname, pos=pos, cigar=cigar, strand=strand)
})
