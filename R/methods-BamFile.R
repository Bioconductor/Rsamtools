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

.normargParam <- function(param, what0)
{
    if (is.null(param))
        param <- ScanBamParam(what=character(0))
    if (!identical(bamFlag(param), scanBamFlag()))
        stop("don't know how to handle a user supplied flag yet")
    flag0 <- scanBamFlag(isUnmappedQuery=FALSE, isDuplicate=FALSE)
    bamFlag(param) <- flag0
    bamWhat(param) <- union(bamWhat(param), what0)
    param
}

### 'bamfile' must be a BamFile object.
### Returns a named list with 1 element per loaded column.
.loadBamCols <- function(bamfile, param, what0)
{
    param <- .normargParam(param, what0)
    res <- scanBam(bamfile, param=param)
    if (length(res) == 0L)
        stop("scanBam() returned a list of length zero")
    ans_names <- names(res[[1L]])
    ans <- lapply(ans_names,
                  function(nm) {
                      tmp <- lapply(unname(res), "[[", nm)
                      tmp1 <- tmp[[1L]]
                      if (is.factor(tmp1)) {
                          ## Fast unlisting of a list of factors that have all
                          ## the same levels in the same order.
                          structure(unlist(tmp),
                                    class="factor", levels=levels(tmp1))
                      } else {
                          do.call(c, tmp)  # doesn't work on list of factors
                      }
                  })
    names(ans) <- ans_names
    ans
}

.loadBamSeqlengths <- function(file, seqlevels)
{
    seqlengths <- scanBamHeader(file)[["targets"]]
    if (is.null(seqlengths))
        return(NULL)
    bad <- setdiff(seqlevels, names(seqlengths))
    if (length(bad) == 0L)
        return(seqlengths)
    bad <- paste(bad, collapse="' '")
    msg <- sprintf("'rname' lengths not in BamFile header; seqlengths not used\n  file: %s\n  missing rname(s): '%s'",
                   path(file), bad)
    warning(msg)
    NULL
}

setMethod(readBamGappedAlignments, "BamFile",
          function(file, index=file, use.names=FALSE, param=NULL)
{
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    what0 <- c("rname", "strand", "pos", "cigar")
    if (use.names)
        what0 <- c(what0, "qname")
    bamcols <- .loadBamCols(file, param, what0)
    seqlengths <- .loadBamSeqlengths(file, levels(bamcols[["rname"]]))
    ans <- GappedAlignments(rname=bamcols$rname, pos=bamcols$pos,
                            cigar=bamcols$cigar, strand=bamcols$strand,
                            names=bamcols$qname, seqlengths=seqlengths)
    if (!is.null(param) && length(bamWhat(param)) != 0L) {
        df <- do.call(DataFrame, bamcols[bamWhat(param)])
        elementMetadata(ans) <- df
    }
    ans
})

setMethod(readBamGappedReads, "BamFile",
          function(file, index=file, use.names=FALSE, param=NULL)
{
    require(ShortRead)  # for the GappedReads() constructor
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    what0 <- c("rname", "strand", "pos", "cigar", "seq")
    if (use.names)
        what0 <- c(what0, "qname")
    bamcols <- .loadBamCols(file, param, what0)
    seqlengths <- .loadBamSeqlengths(file, levels(bamcols[["rname"]]))
    ans <- GappedReads(rname=bamcols$rname, pos=bamcols$pos,
                       cigar=bamcols$cigar, strand=bamcols$strand,
                       qseq=bamcols$seq,
                       names=bamcols$qname, seqlengths=seqlengths)
    if (!is.null(param) && length(bamWhat(param)) != 0L) {
        df <- do.call(DataFrame, bamcols[bamWhat(param)])
        elementMetadata(ans) <- df
    }
    ans
})

