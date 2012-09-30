setMethod(isOpen, "BamFile",
    function(con, rw="")
{
    if (!missing(rw) && rw != "read")
        stop("'rw' must be 'read'")
    .Call(.bamfile_isopen, .extptr(con))
})

BamFile <-
    function(file, index=file, ..., yieldSize=NA_integer_)
{
    .RsamtoolsFile(.BamFile, .normalizePath(file),
                   .normalizePath(index), yieldSize=yieldSize,
                   ...)
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
    if (!isOpen(files)) {
        open(files)
        on.exit(close(files))
    }
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
    h <- scanBamHeader(x)[["targets"]]
    Seqinfo(names(h), unname(h))
})

setMethod(scanBam, "BamFile",
          function(file, index=file, ...,
                   param=ScanBamParam(what=scanBamWhat()))
{
    if (!isOpen(file)) {
        open(file)
        on.exit(close(file))
    }
    if (!missing(index))
        warning("'index' ignored for scanBam,BamFile-method")
    if (!is(param, "ScanBamParam")) {
        msg <- sprintf("'%s' must be a '%s'; was '%s'",
                       "param", "ScanBamParam", class(param))
        stop(msg)
    }
    if (0L != length(bamWhich(param)) && !is.na(yieldSize(file))) {
        msg <- sprintf("'%s' must be '%s' when '%s'",
                       "yieldSize(file)", "NA_integer_",
                       "0 != length(bamWhich(param))")
        stop(msg)
    }
    reverseComplement <- bamReverseComplement(param)
    tmpl <- .scanBam_template(param)
    x <- .io_bam(.scan_bamfile, file, reverseComplement,
                 yieldSize(file), tmpl, param=param)
    .scanBam_postprocess(x, param)
})

setMethod(countBam, "BamFile",
          function (file, index=file, ..., param = ScanBamParam())
{
    if (!isOpen(file)) {
        open(file)
        on.exit(close(file))
    }
    if (!missing(index))
        warning("'index' ignored for countBam,BamFile-method")
    x <- .io_bam(.count_bamfile, file, param=param)
    .countBam_postprocess(x, file, param)
})

setMethod(filterBam, "BamFile",
          function (file, destination, index=file, ...,
                    indexDestination=TRUE,
                    param=ScanBamParam(what=scanBamWhat()))
{
    if (!isOpen(file)) {
        open(file)
        on.exit(close(file))
    }
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
    callGeneric(path(file), destination, ...,
                byQname=byQname, maxMemory=maxMemory)
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "readBamGappedAlignments", "readBamGappedReads" and
### "readBamGappedAlignmentPairs" methods.
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

.normargParam <- function(param, flag0, what0)
{
    if (is.null(param))
        param <- ScanBamParam()
    bamFlag(param) <-
        .combineBamFlagFilters(bamFlag(param, asInteger=TRUE), flag0)
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
    flag0 <- scanBamFlag(isUnmappedQuery=FALSE)
    param <- .normargParam(param, flag0, what0)
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
                                         tag_elt <- tag[[nm]]
                                         ## Fill empty tag with NAs.
                                         if (is.null(tag_elt)) {
                                             count <- length(res_elt[["rname"]])
                                             tag_elt <- rep.int(NA, count)
                                         }
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
        mcols(x) <- df
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
    ans <- GappedAlignments(seqnames=bamcols$rname, pos=bamcols$pos,
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
    ans <- GappedReads(seqnames=bamcols$rname, pos=bamcols$pos,
                       cigar=bamcols$cigar, strand=bamcols$strand,
                       qseq=bamcols$seq, seqlengths=seqlengths)
    .bindExtraData(ans, use.names, param, bamcols)
})

setMethod(readBamGappedAlignmentPairs, "BamFile",
          function(file, index=file, use.names=FALSE, param=NULL)
{
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    if (!is.na(yieldSize(file))) {
        warning("'yieldSize' set to 'NA'")
        yieldSize(file) <- NA_integer_
    }
    flag0 <- scanBamFlag(isPaired=TRUE, hasUnmappedMate=FALSE)
    what0 <- c("flag", "mrnm", "mpos")
    param2 <- .normargParam(param, flag0, what0)
    galn <- readBamGappedAlignments(file, use.names=TRUE, param=param2)
    if (is.null(param)) {
        use.mcols <- FALSE
    } else {
        use.mcols <- c(bamWhat(param), bamTag(param))
    }
    makeGappedAlignmentPairs(galn, use.names=use.names, use.mcols=use.mcols)
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "summarizeOverlaps" methods.
###

.processBamFiles <-
    function(features, reads, mode, ignore.strand, ..., singleEnd=TRUE, param)
{
    mode <- match.fun(mode)

    if ("package:parallel" %in% search() & .Platform$OS.type !=
        "windows" )
        lapply <- parallel::mclapply

    lst <- lapply(seq_along(reads), function(i, reads) {
        bf <- reads[[i]]
        if (!isOpen(bf)) {
            open(bf)
            on.exit(close(bf))
        }
        cnt <- integer(length(features))
        if (singleEnd) {
            while (length(x <- readBamGappedAlignments(bf, param=param))) {
                cnt <- cnt + GenomicRanges:::.dispatchOverlaps(x,
                    features, mode=mode, ignore.strand=ignore.strand)
            }
        } else {
            while (length(x <- readBamGappedAlignmentPairs(bf, param=param))) {
                x <- grglist(x)
                cnt <- cnt + GenomicRanges:::.dispatchOverlaps(x,
                    features, mode=mode, ignore.strand=ignore.strand)
            }
        }
        cnt
    }, reads)

    counts <- do.call(cbind, lst)
    colData <- DataFrame(fileName = reads)
    if (!is.null(names(reads))) {
        rownames(colData) <- names(reads)
    } else {
        rownames(colData) <- sub(".bai$", "", basename(reads))
    }
    SummarizedExperiment(assays=SimpleList(counts=as.matrix(counts)),
                         rowData=features, colData=colData)
}

setMethod("summarizeOverlaps", c("GRanges", "BamFileList"),
    function(features, reads, mode, ignore.strand=FALSE, ..., 
             singleEnd=TRUE, param=ScanBamParam())
{
    .processBamFiles(features, reads, mode, ignore.strand, 
        ..., singleEnd=singleEnd, param=param)
})

setMethod("summarizeOverlaps", c("GRangesList", "BamFileList"),
    function(features, reads, mode, ignore.strand=FALSE, ..., 
             singleEnd=TRUE, param=ScanBamParam())
{
    .processBamFiles(features, reads, mode, ignore.strand, 
        ..., singleEnd=singleEnd, param=param)
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "findSpliceOverlaps" methods.
###

setMethod("findSpliceOverlaps", c("character", "ANY"),
          function(query, subject, ignore.strand=FALSE, ...,
                   param=ScanBamParam(), pairedEnd=FALSE)
{
    findSpliceOverlaps(BamFile(query), subject, ignore.strand, ...,
                       param=param, pairedEnd=pairedEnd)
})

setMethod("findSpliceOverlaps", c("BamFile", "ANY"),
    function(query, subject, ignore.strand=FALSE, ...,
             param=ScanBamParam(), pairedEnd=FALSE)
{
    findSpliceOverlaps(.readRanges(query, param, pairedEnd), subject,
                       ignore.strand, ...)
})

.readRanges <- function(bam, param, pairedEnd)
{
    if (!"XS" %in% bamTag(param))
        bamTag(param) <- c(bamTag(param), "XS")
    if (!pairedEnd)
        reads <- readBamGappedAlignments(bam, param=param)
    else {
        reads <- readGappedAlignmentPairs(path(bam), param=param)
        first_xs <- mcols(first(reads))$XS
        last_xs <- mcols(last(reads))$XS
        if (!is.null(first_xs) && !is.null(last_xs)) {
            xs <- first_xs
            xs[is.na(xs)] <- last_xs[is.na(xs)]
            mcols(reads)$XS <- xs
        }
    }
 
    metadata(reads)$bamfile <- bam
    
    reads
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "coverage" methods.
###

setMethod("coverage", "BamFile",
          function(x, shift=0L, width=NULL, weight=1L, ...)
          coverage(readBamGappedAlignments(x),
                   shift=shift, width=width, weight=weight, ...)
          )
