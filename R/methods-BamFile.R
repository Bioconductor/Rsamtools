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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "readBamGappedAlignments" and "readBamGappedReads" methods.
###

### A "flag filter" is represented as a 'flag' vector of length 2 with names
### keep0 and keep1. The .combineBamFlagFilters() function performs a logical
### AND between 2 "flag filters". It returns a "flag filter".
.combineBamFlagFilters <- function(flagfilterA, flagfilterB)
{
    if (!identical(names(flagfilterA), c("keep0", "keep1"))
     || !identical(names(flagfilterB), c("keep0", "keep1")))
        stop("input must be BAM flag filters")
    ans <- bamFlagAND(flagfilterA, flagfilterB)
    if (!all(bamFlagAsBitMatrix(ans[["keep0"]]) |
             bamFlagAsBitMatrix(ans[["keep1"]])))
        stop("BAM flag filters to combine are incompatible")
    ans
}

.normargParam <- function(param, what0)
{
    if (is.null(param))
        param <- ScanBamParam(what=character(0))
    flag0 <- scanBamFlag(isUnmappedQuery=FALSE, isDuplicate=FALSE)
    bamFlag(param) <- .combineBamFlagFilters(bamFlag(param), flag0)
    bamWhat(param) <- union(bamWhat(param), what0)
    param
}

### 'x' must be a list of length >= 1 where all the elements have the same
### type so combining them is straightforward.
.quickUnlist <- function(x)
{
    x1 <- x[[1L]]
    if (is.factor(x1)) {
        ## Fast unlisting of a list of factors that have all
        ## the same levels in the same order.
        structure(unlist(x), class="factor", levels=levels(x1))
    } else {
        do.call(c, x)  # doesn't work on list of factors
    }
}

### 'bamfile' must be a BamFile object.
### Returns a named list with 1 element per loaded column.
.loadBamCols <- function(bamfile, param, what0)
{
    param <- .normargParam(param, what0)
    res <- unname(scanBam(bamfile, param=param))
    if (length(res) == 0L)
        stop("scanBam() returned a list of length zero")
    ## Extract the "what" cols.
    ans1 <- lapply(bamWhat(param),
                   function(nm) {
                       tmp <- lapply(res, "[[", nm)
                       .quickUnlist(tmp)
                   })
    names(ans1) <- bamWhat(param)
    ## Extract the "tag" cols.
    ans2 <- lapply(bamTag(param),
                   function(nm) {
                       tmp <- lapply(res,
                                     function(res_elt) {
                                         tag <- res_elt[["tag"]]
                                         count <- attr(tag, "count")
                                         tag_elt <- tag[[nm]]
                                         ## Fill empty tag with NAs.
                                         if (is.null(tag_elt))
                                             tag_elt <- rep.int(NA, count)
                                         tag_elt
                                     })
                       .quickUnlist(tmp)
                   })
    names(ans2) <- bamTag(param)
    c(ans1, ans2)
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

### 'x' must be a GappedAlignments object.
.bindExtraData <- function(x, use.names, param, bamcols)
{
    if (use.names)
        names(x) <- bamcols$qname
    if (is.null(param))
        return(x)
    colnames <- c(bamWhat(param), bamTag(param))
    if (length(colnames) != 0L) {
        df <- do.call(DataFrame, bamcols[colnames])
        ## The DataFrame() constructor will mangle the duplicated colnames
        ## to make them unique. We don't want this so we fix them. Note
        ## that we cannot use "colnames<-" for this because it's broken,
        ## but "names<-" seems to work as expected.
        names(df) <- colnames
        elementMetadata(x) <- df
    }
    x
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
                            seqlengths=seqlengths)
    .bindExtraData(ans, use.names, param, bamcols)
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
                       qseq=bamcols$seq, seqlengths=seqlengths)
    .bindExtraData(ans, use.names, param, bamcols)
})


## summarizeOverlaps methods
.processBamFiles <- function(reads, features, mode, ignore.strand, ..., param){
    mode <- match.fun(mode)
    if("package:parallel" %in% search() & .Platform$OS.type != "windows" ){
      lapply <- mclapply}
    reads <- path(reads)
    lst <- lapply(reads,
                  function(bf) {
                    x <- readGappedAlignments(bf, param=param)
                    GenomicRanges:::.dispatch(x, features, mode=mode, 
                              ignore.strand=ignore.strand)
                  })
    counts <- do.call(cbind, lst)
    colData <- DataFrame(fileName = reads)
    rownames(colData) <- sub(".bai$", "", basename(reads))
    SummarizedExperiment(assays=SimpleList(counts=as.matrix(counts)),
                         rowData=features, colData=colData)
}

setMethod("summarizeOverlaps", c("BamFileList", "GRanges"),
    function(reads, features, 
             mode, 
             ignore.strand = FALSE, ..., param = ScanBamParam())
{
    .processBamFiles(reads, features, mode, ignore.strand, ..., param=param)
})

setMethod("summarizeOverlaps", c("BamFileList", "GRangesList"),
    function(reads, features, 
             mode, 
             ignore.strand = FALSE, ..., param = ScanBamParam())
{
    .processBamFiles(reads, features, mode, ignore.strand, ..., param=param)
})


