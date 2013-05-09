setMethod(isOpen, "BamFile",
    function(con, rw="")
{
    if (!missing(rw) && rw != "read")
        stop("'rw' must be 'read'")
    .Call(.bamfile_isopen, .extptr(con))
})

setMethod(isIncomplete, "BamFile",
    function(con)
{
    .Call(.bamfile_isincomplete, .extptr(con))
})

BamFile <-
    function(file, index=file, ..., yieldSize=NA_integer_, 
             obeyQname=FALSE)
{
    .RsamtoolsFile(.BamFile, .normalizePath(file),
                   .normalizePath(index), yieldSize=yieldSize,
                   obeyQname=obeyQname, ...)
}

open.BamFile <-
    function(con, ...)
{
    tryCatch({
        .io_check_exists(path(con))
        index <- sub("\\.bai$", "", index(con))
        con$.extptr <- .Call(.bamfile_open, path(con), index, "rb")
    }, error=function(err) {
        stop("failed to open BamFile: ", conditionMessage(err))
    })
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

setMethod(obeyQname, "BamFile",
    function(object, ...)
{
    object$obeyQname
})

setReplaceMethod("obeyQname", "BamFile", 
    function(object, ..., value)
{
    if (1L != length(value))
        stop("'value' must be length 1")
    object$obeyQname <- value
    object
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
                 yieldSize(file), tmpl, obeyQname(file), param=param)
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
                    filter=FilterRules(),
                    indexDestination=TRUE,
                    param=ScanBamParam(what=scanBamWhat()))
{
    if (!isOpen(file)) {
        open(file)
        on.exit(close(file))
    }
    if (missing(destination))
        stop(sprintf("'%s' missing with no default", "destination"))
    if (!is(param, "ScanBamParam"))
        stop(sprintf("'%s' must be a '%s'; was '%s'",
                     "param", "ScanBamParam", class(param)))

    seqlengths <- seqlengths(file)
    yieldSize <- yieldSize(file)
    if (is.na(yieldSize))
        yieldSize <- 1000000L
    param <- .filterBam_preprocess(file, param)
    destination <- .normalizePath(destination)
    dest <- .Call(.bamfile_open, destination, path(file), "wb")
    n_tot <- 0L

    repeat {
        buf <- .io_bam(.prefilter_bamfile, file, param=param, yieldSize,
                       obeyQname(file))
        if (0L == .Call(.bambuffer_length, buf))
            break;

        if (length(filter)) {
            tmpl <- .scanBam_template(param)
            reverseComplement <- bamReverseComplement(param)
            ans <- .io_bam(.bambuffer_parse, file, param=param,
                           buf, reverseComplement, tmpl)
            ans <- .loadBamColsFromScan(unname(ans), param)
            ans <- DataFrame(ans)
            ans <- eval(filter, ans)
        } else ans <- TRUE

        n_tot <- n_tot + .Call(.bambuffer_write, buf, dest, ans)
    }

    .Call(.bamfile_close, dest)
    if (indexDestination)
        indexBam(destination)
    destination
})

setMethod(indexBam, "BamFile", function(files, ...) {
    indexBam(path(files), ...)
})

setMethod(sortBam, "BamFile",
    function(file, destination, ..., byQname=FALSE, maxMemory=512)
{
    sortBam(path(file), destination, ...,
                byQname=byQname, maxMemory=maxMemory)
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "readGAlignmentsFromBam", "readGappedReadsFromBam",
### "readGAlignmentPairsFromBam", and "readGAlignmentsListFromBam" methods.
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
.loadBamColsFromScan <- function(res, param)
{
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

.loadBamCols <- function(bamfile, param, what0, ...)
{
    flag0 <- scanBamFlag(isUnmappedQuery=FALSE)
    param <- .normargParam(param, flag0, what0)
    res <- unname(scanBam(bamfile, ..., param=param))
    if (length(res) == 0L)
        stop("scanBam() returned a list of length zero")
    .loadBamColsFromScan(res, param)
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

### 'x' must be a GAlignments object.
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

setMethod(readGAlignmentsFromBam, "BamFile",
          function(file, index=file, ..., use.names=FALSE, param=NULL)
{
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    what0 <- c("rname", "strand", "pos", "cigar")
    if (use.names)
        what0 <- c(what0, "qname")
    bamcols <- .loadBamCols(file, param, what0, ...)
    seqlengths <- .loadBamSeqlengths(file, levels(bamcols[["rname"]]))
    ans <- GAlignments(seqnames=bamcols$rname, pos=bamcols$pos,
                       cigar=bamcols$cigar, strand=bamcols$strand,
                       seqlengths=seqlengths)
    .bindExtraData(ans, use.names, param, bamcols)
})

setMethod(readGappedReadsFromBam, "BamFile",
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

setMethod(readGAlignmentPairsFromBam, "BamFile",
          function(file, index=file, use.names=FALSE, param=NULL)
{
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    if (!isTRUE(obeyQname(file))  && !is.na(yieldSize(file))) {
        warning("'yieldSize' set to 'NA'")
        yieldSize(file) <- NA_integer_
    }
    flag0 <- scanBamFlag(isPaired=TRUE, hasUnmappedMate=FALSE)
    what0 <- c("flag", "mrnm", "mpos")
    param2 <- .normargParam(param, flag0, what0)
    galn <- readGAlignmentsFromBam(file, use.names=TRUE, param=param2)
    if (is.null(param)) {
        use.mcols <- FALSE
    } else {
        use.mcols <- c(bamWhat(param), bamTag(param))
    }
    makeGAlignmentPairs(galn, use.names=use.names, use.mcols=use.mcols)
})

setMethod(readGAlignmentsListFromBam, "BamFile", 
    function(file, index=file, ..., use.names=FALSE, param=ScanBamParam(),
             group.as.pairs=TRUE)
{
    if (yieldSize(file) && !obeyQname(file))
        stop("'obeyQname(file)' must be TRUE when 'yieldSize(file)' != NA")
    if (group.as.pairs) {
        use.mcols=c(bamWhat(param), bamTag(param))
        bamWhat(param) <- union(bamWhat(param), c("flag", "mrnm", "mpos"))
        gal <- readGAlignmentsFromBam(file, use.names=TRUE, param=param)
        groupAsPairs(gal, use.mcols=use.mcols) 
    } else {
        gal <- readGAlignmentsFromBam(file, use.names=TRUE, param=param)
        splitAsList(gal, factor(names(gal)))  ## may not be ordered
    }
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "groupAsPairs" methods.
###

setMethod(groupAsPairs, "GAlignmentsList", 
    function(x, ...) 
{
    callGeneric(unlist(x, use.names=FALSE), ...) 
})

setMethod(groupAsPairs, "GAlignments", 
    function(x, ..., use.mcols=colnames(mcols(x))) 
{
    if (length(use.mcols) > 0L) {
        if (!all(use.mcols %in% colnames(mcols(x))))
            stop("'use.mcols' must be a subset of 'colnames(mcols(x))'")
    }
    .groupAsPairs(x, use.mcols=use.mcols)
})

.groupAsPairs <- function(x, use.mcols)
{
    ## Potential mates.
    mate <- findMateAlignment(x)
    x_is_first <- .isFirstSegment.GAlignments(x)
    x_is_last <- .isLastSegment.GAlignments(x)
    first_idx <- which(!is.na(mate) & x_is_first)
    last_idx <- mate[first_idx]
    .checkMates(mate, x_is_first, x_is_last, first_idx, last_idx)

    ## Check the 0x2 bit (isProperPair).
    x_is_proper <- as.logical(bamFlagAsBitMatrix(mcols(x)$flag,
                              bitnames="isProperPair"))
    ans_is_proper <- x_is_proper[first_idx]

    ## Check pairs for discordant seqnames or strand.
    x_is_discordant <- (as.character(seqnames(x)[first_idx]) !=
                        as.character(seqnames(x)[last_idx])) | 
                       (as.character(strand(x)[first_idx]) ==
                        as.character(strand(x)[last_idx]))
    keep <- which(!x_is_discordant)
    first_idx <- first_idx[keep]
    last_idx <- last_idx[keep]
    ## FIXME: is this ever FALSE?
    ans_is_proper <- ans_is_proper[keep]

    ## Assemble GAlignmentsList.
    pairs_idx <- c(rbind(first_idx, last_idx))
    other_idx <- seq_along(x)[-pairs_idx]
    mcols(x)$paired <- seq_along(x) %in% pairs_idx
    mcols(x) <- mcols(x)[c(use.mcols, "paired")] 
    other_nms <- names(x)[other_idx]
    other_elt <- splitAsList(other_nms, factor(other_nms))
    widths <- c(rep(2, length(first_idx)), elementLengths(other_elt))
    relist(x[c(pairs_idx, other_idx)], PartitioningByEnd(cumsum(widths)))
}

.checkMates <- function(mate, x_is_first, x_is_last, first_idx, last_idx)
{
    ## Fundamental property of the 'mate' vector: it's a permutation of order
    ## 2 and with no fixed point on the set of indices for which 'mate' is
    ## not NA.
    ## Check there are no fixed points.
    if (!all(first_idx != last_idx))
        stop("findMateAlignment() returned an invalid 'mate' vector")
    ## Check order 2 (i.e. permuting a 2nd time brings back the original
    ## set of indices).
    if (!identical(mate[last_idx], first_idx))
        stop("findMateAlignment() returned an invalid 'mate' vector")
    ## One more sanity check.
    if (!all(x_is_last[last_idx]))
        stop("findMateAlignment() returned an invalid 'mate' vector")
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "summarizeOverlaps" methods.
###

setMethod("summarizeOverlaps", c("GRanges", "BamFile"),
    function(features, reads, mode, ignore.strand=FALSE, ..., 
             inter.feature=TRUE, singleEnd=TRUE, fragments=TRUE,
             param=ScanBamParam())
{
    .dispatchBamFiles(features, BamFileList(reads), mode, ignore.strand, ..., 
                      inter.feature=inter.feature, singleEnd=singleEnd, 
                      fragments=fragments, param=param)
})

setMethod("summarizeOverlaps", c("GRangesList", "BamFile"),
    function(features, reads, mode, ignore.strand=FALSE, ..., 
             inter.feature=TRUE, singleEnd=TRUE, fragments=TRUE, 
             param=ScanBamParam())
{
    .dispatchBamFiles(features, BamFileList(reads), mode, ignore.strand, ..., 
                      inter.feature=inter.feature, singleEnd=singleEnd, 
                      fragments=fragments, param=param)
})

.dispatchOverlaps <- GenomicRanges:::.dispatchOverlaps
.countWithYieldSize <- function(FUN, features, bf, mode, ignore.strand, 
                                inter.feature, param)
{
    if (is.na(yieldSize(bf))) {
        x <- FUN(bf, param=param)
        .dispatchOverlaps(features, x, mode, ignore.strand, inter.feature)
    } else {
        if (!isOpen(bf)) {
            open(bf)
            on.exit(close(bf))
        }
        ct <- integer(length(features))
        while (length(x <- FUN(bf, param=param))) {
            ct <- ct + .dispatchOverlaps(features, x, mode, ignore.strand,
                                         inter.feature) 
        }
        ct
    }
}

.dispatchBamFiles <-
    function(features, reads, mode, ignore.strand, ..., 
             inter.feature=TRUE, singleEnd=TRUE, fragments=TRUE,
             param=ScanBamParam())
{
    if ("package:parallel" %in% search() & .Platform$OS.type != "windows")
        lapply <- parallel::mclapply

    cts <- lapply(seq_along(reads), 
               function(i, reads, mode) {
                   bf <- reads[[i]]
                   if (singleEnd) {
                       ## single-end
                       FUN <- readGAlignmentsFromBam
                   } else {
                       if (fragments) { 
                           if (isTRUE(obeyQname(bf))) 
                               ## paired-end, sorted by qname
                               FUN <- readGAlignmentsListFromBam
                           else
                               stop("when 'fragments=TRUE' Bam files must be ",
                                    "sorted by qname ('obeyQname=TRUE')")
                       } else {
                           ## paired-end, not sorted by qname
                           FUN <- readGAlignmentPairsFromBam
                       }
                   }
                   .countWithYieldSize(FUN, features, bf, mode, ignore.strand, 
                                       inter.feature=inter.feature, param=param) 
               }, reads, mode=match.fun(mode))

    counts <- as.matrix(do.call(cbind, cts))
    colData <- DataFrame(fileName=reads)
    if (!is.null(names(reads)))
        rownames(colData) <- names(reads)
    else
        rownames(colData) <- sub(".bai$", "", basename(reads))
    SummarizedExperiment(assays=SimpleList(counts=counts),
                         rowData=features, colData=colData)
}

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
        reads <- readGAlignmentsFromBam(bam, param=param)
    else {
        reads <- readGAlignmentPairsFromBam(path(bam), param=param)
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
          function(x, shift=0L, width=NULL, weight=1L, ...,
                   param = ScanBamParam())
          coverage(readGAlignmentsFromBam(x, param = param),
                   shift=shift, width=width, weight=weight, ...)
          )

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "quickCountBam" methods.
###

setMethod(quickCountBam, "BamFile",
    function(file, ..., param=ScanBamParam(), mainGroupsOnly=FALSE)
{
    what0 <- c("qname", "flag")
    if (length(bamWhat(param)) != 0L)
        warning("bamWhat component of supplied 'param' was ignored")
    bamWhat(param) <- what0

    res <- scanBam(file, param=param)
    res0 <- res[[1L]]
    if (length(res) != 1L) {
        res0[["qname"]] <- do.call(c, lapply(res, "[[", "qname"))
        res0[["flag"]] <- do.call(c, lapply(res, "[[", "flag"))
    } 

    quickCountBam(res0, param=param, mainGroupsOnly=mainGroupsOnly)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "show" methods
###

setMethod(show, "BamFile", function(object) {
    callNextMethod()
    cat("obeyQname:", obeyQname(object), "\n")
})
