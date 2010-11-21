.BamFile <- setRefClass("BamFile",
    fields=list(
      .extptr="externalptr",
      bamPath="character",
      bamIndex="character"))

setGeneric("isOpen")

setGeneric("openBam",
           function(bamPath, bamIndex=bamPath, ...)
           standardGeneric("openBam"), signature="bamPath")

.extptr <- function(object) object$.extptr

bamPath <- function(object) object$bamPath

bamIndex <- function(object) object$bamIndex

setMethod(isOpen, "BamFile",
    function(con, rw="")
{
    if (!missing(rw) && rw == "read")
        stop("'rw' must be 'read'")
    .Call(.bamfile_isopen, con$.extptr)
})

.openBam <-
    function(bamPath, bamIndex=bamPath)
{
    .io_bam_check_exists(bamPath)
    bamPath <- .normalizePath(bamPath)
    bamIndex <- .normalizePath(bamIndex)
    extptr <- .Call(.bamfile_open, bamPath, bamIndex, "rb")
    .BamFile$new(.extptr=extptr, bamPath=bamPath,
                 bamIndex=bamIndex)
}

setMethod(openBam, "character",
          function(bamPath, bamIndex=bamPath, ...)
{
    .openBam(bamPath, bamIndex)
})

setMethod(openBam, "BamFile",
          function(bamPath, bamIndex=bamPath, ...)
{
    bf <- bamPath                      # convenience
    callGeneric(bamPath(bf), bamIndex(bf), ...)
})

close.BamFile <-
    function(con, ...)
{
    if (!isOpen(con))
        stop("isOpen(<BamFile>) is not 'TRUE'")
    con$.extptr <- .Call(.bamfile_close, .extptr(con))
    invisible(con)
}

setMethod(show, "BamFile", function(object) {
    .ppath <- function(tag, filepath)
    {
        wd <- options('width')[[1]] - nchar(tag) - 6
        if (0L == length(filepath) || nchar(filepath) < wd)
            return(sprintf("%s: %s\n", tag, filepath))
        bname <- basename(filepath)
        wd1 <- wd - nchar(bname)
        dname <- substr(dirname(filepath), 1, wd1)
        sprintf("%s: %s...%s%s\n",
                tag, dname, .Platform$file.sep, bname)
    }
    cat("class:", class(object), "\n")
    cat(.ppath("bamPath", bamPath(object)))
    cat(.ppath("bamIndex", bamIndex(object)))
    cat("isOpen:", isOpen(object), "\n")
})

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
                 bamPath(file), bamIndex(file), "rb",
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
    callGeneric(bamPath(files), ...)
})

setMethod(sortBam, "BamFile",
    function(file, destination, ..., byQname=FALSE, maxMemory=512)
{
    callGeneric(bamPath(files), destination, ...,
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
