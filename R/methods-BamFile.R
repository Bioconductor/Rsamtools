setMethod(isOpen, "BamFile",
    function(con, rw="")
{
    if (!missing(rw) && rw != "read")
        stop("'rw' must be 'read'")
    .Call(.bamfile_isopen, .extptr(con))
})

BamFile <-
    function(file, index=file)
{
    .RsamtoolsFile(.BamFile, .normalizePath(file),
                   .normalizePath(index))
}

open.BamFile <-
    function(con, ...)
{
    .io_check_exists(path(con))
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
    if (!isOpen(files))
        stop("BamFile is not open")
    header <- .Call(.read_bamfile_header, .extptr(files))
    text <- strsplit(header[["text"]], "\n")[[1]]
    tag <- sub("^(@[A-Z]{2}).*", "\\1", text)
    text <- strsplit(sub("^@[A-Z]{2}\t(.*)", "\\1", text), "\t")
    names(text) <- tag
    header[["text"]] <- text
    header
})

setMethod(seqinfo, "BamFile",
          function(x)
{
    if (!isOpen(x))
        stop("BamFile is not open")
    h <- scanBamHeader(x)[["targets"]]
    Seqinfo(names(h), unname(h))
})

setMethod(scanBam, "BamFile",
          function(file, index=file, ..., param=ScanBamParam())
{
    if (!isOpen(file))
        stop("BamFile not open")
    if (!missing(index))
        warning("'index' ignored for scanBam,BamFile-method")
    if (!is(param, "ScanBamParam")) {
        msg <- sprintf("'%s' must be a '%s'; was '%s'",
                       "param", "ScanBamParam", class(param))
        stop(msg)
    }
    reverseComplement <- bamReverseComplement(param)
    tmpl <- .scanBam_template(param)
    x <- .io_bam(.scan_bamfile, file, param=param, reverseComplement,
                 tmpl)
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

### 'bamfile' must be a BamFile object.
### Returns either NULL or a list which names are guaranted to be 'what'.
.loadBamCols <- function(bamfile, what, which)
{
    param <- ScanBamParam(
                 flag=scanBamFlag(isUnmappedQuery=FALSE, isDuplicate=FALSE),
                 what=what,
                 which=which)
    bam <- scanBam(bamfile, param=param)
    ## unlist(list(factor())) returns integer(0), so exit early if all
    ## values are empty
    bam <- bam[sapply(bam, function(x) length(x$rname) != 0L)]
    if (length(bam) == 0L)
        return(NULL)
    if (length(bam) == 1L)
        return(bam[[1L]][what])
    ans <- lapply(what,
                  function(field) {
                      tmp <- unname(lapply(bam, "[[", field))
                      if (is.factor(tmp[[1L]]))
                          unlist(tmp)  # doesn't work on list of XStringSets
                      else
                          do.call(c, tmp)  # doesn't work on list of factors
                  })
    names(ans) <- what
    ans
}

setMethod(readBamGappedAlignments, "BamFile",
          function(file, index=file, use.names=FALSE, ..., which)
{
    if (missing(which))
        which <- RangesList()
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    what <- c("rname", "strand", "pos", "cigar")
    if (use.names)
        what <- c(what, "qname")
    seqlengths <- scanBamHeader(file)[["targets"]]
    bamcols <- .loadBamCols(file, what, which)
    if (is.null(bamcols)) {
        ans <- GappedAlignments(rname=factor(levels=names(seqlengths)),
                                seqlengths=seqlengths)
        return(ans)
    }
    ## 'bamcols' names must be valid GappedAlignments() arg names.
    what[match("qname", what)] <- "names"
    names(bamcols) <- what
    seqlevels <- levels(bamcols[["rname"]])
    if (all(seqlevels %in% names(seqlengths))) {
        bamcols$seqlengths <- seqlengths
    } else {
        bad <- paste(seqlevels[!(seqlevels %in% names(seqlengths))],
                     collapse="' '")
        msg <- sprintf("'rname' lengths not in BamFile header; seqlengths not used\n  file: %s\n  missing rname(s): '%s'",
                       path(file), bad)
        warning(msg)
    }
    do.call(GappedAlignments, bamcols)
})

setMethod(readBamGappedReads, "BamFile",
          function(file, index=file, use.names=FALSE, ..., which)
{
    require(ShortRead)  # for the GappedReads() constructor
    if (missing(which))
        which <- RangesList()
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    what <- c("rname", "strand", "pos", "cigar", "seq")
    if (use.names)
        what <- c(what, "qname")
    seqlengths <- scanBamHeader(file)[["targets"]]
    bamcols <- .loadBamCols(file, what, which)
    if (is.null(bamcols)) {
        ans <- GappedReads(rname=factor(levels=names(seqlengths)),
                           seqlengths=seqlengths)
        return(ans)
    }
    ## 'bamcols' names must be valid GappedReads() arg names.
    what[match("qname", what)] <- "names"
    what[match("seq", what)] <- "qseq"
    names(bamcols) <- what
    seqlevels <- levels(bamcols[["rname"]])
    if (all(seqlevels %in% names(seqlengths))) {
        bamcols$seqlengths <- seqlengths
    } else {
        bad <- paste(seqlevels[!seqlevels %in% names(seqlengths)],
                     collapse="' '")
        msg <- sprintf("'rname' lengths not in BamFile header; seqlengths not used\n  file: %s\n  missing rname(s): '%s'",
                       path(file), bad)
        warning(msg)
    }
    do.call(GappedReads, bamcols)
})

