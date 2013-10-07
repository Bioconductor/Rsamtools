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
             obeyQname=FALSE, asMates=FALSE)
{
    .RsamtoolsFile(.BamFile, .normalizePath(file),
                   .normalizePath(index), yieldSize=yieldSize,
                   obeyQname=obeyQname, asMates=asMates, ...)
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

setMethod(asMates, "BamFile",
    function(object, ...)
{
    object$asMates
})

setReplaceMethod("asMates", "BamFile", 
    function(object, ..., value)
{
    if (1L != length(value))
        stop("'value' must be length 1")
    object$asMates <- value
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
    if (obeyQname(file) && asMates(file))
        warning("'obeyQname=TRUE' ignored when 'asMates=TRUE'")
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
    if (!asMates(file))
        bamWhat(param) <- setdiff(bamWhat(param), c("partition", "mates")) 
    reverseComplement <- bamReverseComplement(param)
    tmpl <- .scanBam_template(param)
    x <- .io_bam(.scan_bamfile, file, reverseComplement,
                 yieldSize(file), tmpl, obeyQname(file), 
                 asMates(file), param=param)
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

.filterBam_FilterRules <-
    function(file, destination, filter, param)
{
    which <- unlist(bamWhich(param))
    nRange <- length(which)
    if (nRange)                         # yield by range
        iRange <- 1L

    yieldSize <- yieldSize(file)
    if (is.na(yieldSize))
        yieldSize <- 1000000L

    tmpl <- .scanBam_template(param)
    reverseComplement <- bamReverseComplement(param)
    partitionAsWidth <- FALSE

    dest <- .Call(.bamfile_open, destination, path(file), "wb")
    n_tot <- 0L
    repeat {
        if (nRange) {                   # by range
            if (iRange > nRange)
                break
            which0 <- IRangesList(which[iRange])
            names(which0) <- names(which)[iRange]
            param <- initialize(param, which=which0)
            iRange <- iRange + 1L
        }
        buf <- .io_bam(.prefilter_bamfile, file, param=param,
                       yieldSize, obeyQname(file), asMates(file))
        if (0L == .Call(.bambuffer_length, buf))
            if (nRange) {
                next
            } else break

        ans <- .io_bam(.bambuffer_parse, file, param=param,
                       buf, reverseComplement, partitionAsWidth, tmpl)
        ans <- DataFrame(.load_bamcols_from_scanBam_res(ans, param))
        ans <- eval(filter, ans)
        n_tot <- n_tot + .Call(.bambuffer_write, buf, dest, ans)
    }

    .Call(.bamfile_close, dest)
    destination
}

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
    param <- .filterBam_preprocess(file, param)
    destination <- .normalizePath(destination)

    if (length(filter))
        .filterBam_FilterRules(file, param=param, destination, filter)
    else
        .io_bam(.filter_bamfile, file, param=param, destination, "wb")

    if (indexDestination) {
        if (asMates(file)) {
            ## FIXME: filtering by mates requires expensive re-sort!
            fl <- tempfile()
            file.rename(destination, fl)
            sortBam(fl, destination)
            file.rename(paste0(destination, ".bam"), destination)
        }
        indexBam(destination)
    }
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
.load_bamcols_from_scanBam_res <- function(res, param, with.which_label=FALSE)
{
    if (!isTRUEorFALSE(with.which_label))
        stop("'with.which_label' must be TRUE or FALSE")
    which_label <- names(res)
    names(res) <- NULL
    ## Compute the "which_label" col.
    ans3 <- list()
    if (with.which_label) {
        if (is.null(which_label)) {
            warning("'which_label' is ignored when 'param' is missing or ",
                    "doesn't have a 'which' component")
        } else {
            ## We want unique labels.
            if (anyDuplicated(which_label))
                which_label <- paste0(which_label, ".", seq_along(which_label))
            ## Currently scanBam() should always output a 'res' of length != 0
            ## but we test anyway just in case one day this changes.
            if (length(res) == 0L) {
                run_lens <- integer(0)
            } else {
                run_lens <- sapply(res, function(x) length(x[[1L]]))
            }
            which_label <- Rle(factor(which_label, levels=which_label),
                               run_lens)
            ans3 <- list(which_label=which_label)
        }
    }
    ## Extract the "what" cols.
    ans1 <- lapply(setNames(bamWhat(param), bamWhat(param)),
                   function(nm) {
                       tmp <- lapply(res, "[[", nm)
                       .quickUnlist(tmp)
                   })
    ## Extract the "tag" cols.
    ans2 <- lapply(setNames(bamTag(param), bamTag(param)),
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
    ## Put all the cols together and return them.
    c(ans1, ans2, ans3)
}

.load_bamcols_from_bamfile <- function(bamfile, param, what0, ...,
                                       with.which_label=FALSE)
{
    flag0 <- scanBamFlag(isUnmappedQuery=FALSE)
    param <- .normargParam(param, flag0, what0)
    res <- scanBam(bamfile, ..., param=param)
    if (length(res) == 0L)  # should never happen
        stop("scanBam() returned a list of length zero")
    .load_bamcols_from_scanBam_res(res, param,
                                   with.which_label=with.which_label)
}

.load_seqlengths_from_bamfile <- function(file, seqlevels)
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
.bindExtraData <- function(x, use.names, param, bamcols,
                           with.which_label=FALSE)
{
    if (use.names)
        names(x) <- bamcols$qname
    if (is.null(param))
        return(x)
    colnames <- c(bamWhat(param), bamTag(param))
    if (with.which_label)
        colnames <- c(colnames, "which_label")
    if (length(colnames) != 0L) {
        df <- do.call(DataFrame, bamcols[colnames])
        ## Sadly, the DataFrame() constructor is mangling the duplicated
        ## colnames to make them unique. Since we of course don't want this,
        ## we need to fix them.
        colnames(df) <- colnames
        mcols(x) <- df
    }
    x
}

setMethod(readGAlignmentsFromBam, "BamFile",
          function(file, index=file, ..., use.names=FALSE, param=NULL,
                   with.which_label=FALSE)
{
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    if (is.null(param))
        param <- ScanBamParam()
    if (!asMates(file))
        bamWhat(param) <- setdiff(bamWhat(param), c("partition", "mates"))
    what0 <- c("rname", "strand", "pos", "cigar")
    if (use.names)
        what0 <- c(what0, "qname")
    bamcols <- .load_bamcols_from_bamfile(file, param, what0, ...,
                                          with.which_label=with.which_label)
    seqlengths <- .load_seqlengths_from_bamfile(file, levels(bamcols[["rname"]]))
    ans <- GAlignments(seqnames=bamcols$rname, pos=bamcols$pos,
                       cigar=bamcols$cigar, strand=bamcols$strand,
                       seqlengths=seqlengths)
    .bindExtraData(ans, use.names, param, bamcols,
                   with.which_label=with.which_label)
})

setMethod(readGappedReadsFromBam, "BamFile",
          function(file, index=file, use.names=FALSE, param=NULL,
                   with.which_label=FALSE)
{
    require(ShortRead)  # for the GappedReads() constructor
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    if (is.null(param))
        param <- ScanBamParam()
    if (!asMates(file))
        bamWhat(param) <- setdiff(bamWhat(param), c("partition", "mates"))
    what0 <- c("rname", "strand", "pos", "cigar", "seq")
    if (use.names)
        what0 <- c(what0, "qname")
    bamcols <- .load_bamcols_from_bamfile(file, param, what0,
                                          with.which_label=with.which_label)
    seqlengths <- .load_seqlengths_from_bamfile(file, levels(bamcols[["rname"]]))
    ans <- GappedReads(seqnames=bamcols$rname, pos=bamcols$pos,
                       cigar=bamcols$cigar, strand=bamcols$strand,
                       qseq=bamcols$seq, seqlengths=seqlengths)
    .bindExtraData(ans, use.names, param, bamcols,
                   with.which_label=with.which_label)
})

setMethod(readGAlignmentPairsFromBam, "BamFile",
          function(file, index=file, use.names=FALSE, param=NULL,
                   with.which_label=FALSE)
{
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    if (is.null(param))
        param <- ScanBamParam()
    if (!asMates(file))
        bamWhat(param) <- setdiff(bamWhat(param), c("partition", "mates"))
    if (!is.na(yieldSize(file))) {
        warning("'yieldSize' set to 'NA'", immediate.=TRUE)
        yieldSize(file) <- NA_integer_
    }
    flag0 <- scanBamFlag(isPaired=TRUE, hasUnmappedMate=FALSE)
    what0 <- c("flag", "mrnm", "mpos")
    param2 <- .normargParam(param, flag0, what0)
    galn <- readGAlignmentsFromBam(file, use.names=TRUE, param=param2,
                                   with.which_label=with.which_label)
    if (is.null(param)) {
        use.mcols <- FALSE
    } else {
        use.mcols <- c(bamWhat(param), bamTag(param))
        if (with.which_label)
            use.mcols <- c(use.mcols, "which_label")
    }
    makeGAlignmentPairs(galn, use.names=use.names, use.mcols=use.mcols)
})

setMethod(readGAlignmentsListFromBam, "BamFile", 
    function(file, index=file, ..., use.names=FALSE, param=ScanBamParam())
{
    ## FIXME: turn on when summarizeOverlaps() is hooked into new pairing algo
    if (!asMates(file)) {
    #    warning("'asMates' should be true; use readGAlignments() for ",
    #            "single-end data.")
    #    bamWhat(param) <- setdiff(bamWhat(param), c("partition", "mates"))
    } else {
        bamWhat(param) <- union("mates", bamWhat(param))
    }
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")

    ## required for GAlignments
    what0 <- c("rname", "strand", "pos", "cigar", "partition")
    if (use.names)
        what0 <- c(what0, "qname")
    .matesFromBam(file, use.names=use.names, param, what0) 
})

.matesFromBam <- function(file, use.names, param, what0)
{
    bamcols <- .load_bamcols_from_bamfile(file, param, what0)
    seqlengths <- .load_seqlengths_from_bamfile(file, levels(bamcols$rname))
    gal <- GAlignments(seqnames=bamcols$rname, pos=bamcols$pos,
                       cigar=bamcols$cigar, strand=bamcols$strand,
                       seqlengths=seqlengths)
    res <- .bindExtraData(gal, use.names=FALSE, param, bamcols)
    if (asMates(file))
        pbw <- PartitioningByWidth(bamcols$partition)
    else
        pbw <- PartitioningByWidth(rep(1, length(gal)))

    galist <- relist(res, pbw)
    ## FIXME: no need to split when 'partition' is fixed
    if (use.names) {
        nms <- unique(relist(bamcols$qname, pbw))
        names(galist) <- nms
    } 
    galist
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "summarizeOverlaps" methods.
###

setMethod("summarizeOverlaps", c("GRanges", "BamFile"),
    function(features, reads, mode, ignore.strand=FALSE, ..., 
             inter.feature=TRUE, singleEnd=TRUE, fragments=TRUE,
             param=ScanBamParam())
{
    .checkArgs(fragments=fragments, singleEnd=singleEnd, bam=reads)
    .dispatchBamFiles(features, BamFileList(reads), mode, ignore.strand, ..., 
                      inter.feature=inter.feature, singleEnd=singleEnd, 
                      fragments=fragments, param=param)
})

setMethod("summarizeOverlaps", c("GRangesList", "BamFile"),
    function(features, reads, mode, ignore.strand=FALSE, ..., 
             inter.feature=TRUE, singleEnd=TRUE, fragments=TRUE, 
             param=ScanBamParam())
{
    .checkArgs(fragments=fragments, singleEnd=singleEnd, bam=reads)
    .dispatchBamFiles(features, BamFileList(reads), mode, ignore.strand, ..., 
                      inter.feature=inter.feature, singleEnd=singleEnd, 
                      fragments=fragments, param=param)
})

.checkArgs <- function(fragments, singleEnd, bam)
{
    if (singleEnd) {
        if (asMates(bam))
            stop("cannot specify both singleEnd=TRUE and asMates=TRUE")
    } else {
        if (fragments && !obeyQname(bam))
            stop("when 'fragments=TRUE' Bam files must be ",
                 "sorted by qname ('obeyQname=TRUE')")
    }
}

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

.getReadFunction <- function(singleEnd, fragments, obeyQname)
{
    if (singleEnd) {
        FUN <- readGAlignmentsFromBam
    } else {
        if (fragments)
            if (obeyQname)
                FUN <- readGAlignmentsListFromBam
            else
                stop("when 'fragments=TRUE' Bam files must be ",
                     "sorted by qname ('obeyQname=TRUE')")
        else
            FUN <- readGAlignmentPairsFromBam
    }
    FUN
}

.dispatchBamFiles <-
    function(features, reads, mode, ignore.strand, ...,
             count.mapped.reads=FALSE,
             inter.feature=TRUE, singleEnd=TRUE, fragments=TRUE,
             param=ScanBamParam())
{
    FUN <- .getReadFunction(singleEnd, fragments,
                            isTRUE(all(obeyQname(reads))))

    if ("package:parallel" %in% search() & .Platform$OS.type != "windows")
        lapply <- parallel::mclapply

    cts <- lapply(setNames(seq_along(reads), names(reads)), 
               function(i, FUN, reads, features, mode, ignore.strand, 
                        inter.feature, param) {
                   bf <- reads[[i]]
                   .countWithYieldSize(FUN, features, bf, mode, ignore.strand, 
                                       inter.feature, param) 
               }, FUN, reads, features, mode=match.fun(mode), ignore.strand, 
               inter.feature, param
           ) 

    counts <- as.matrix(do.call(cbind, cts))
    if (count.mapped.reads) {
        countBam <- countBam(reads)
        flag <- scanBamFlag(isUnmappedQuery=FALSE)
        param <- ScanBamParam(flag=flag, what="seq")
        colData <- DataFrame(countBam[c("records", "nucleotides")],
                             mapped=countBam(reads, param=param)$records,
                             row.names=colnames(counts))
    } else {
        colData <- DataFrame(row.names=colnames(counts))
    }
    SummarizedExperiment(assays=SimpleList(counts=counts),
                         rowData=features, colData=colData)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "findSpliceOverlaps" methods.
###

setMethod("findSpliceOverlaps", c("character", "ANY"),
          function(query, subject, ignore.strand=FALSE, ...,
                   param=ScanBamParam(), singleEnd=TRUE)
{
    findSpliceOverlaps(BamFile(query), subject, ignore.strand, ...,
                       param=param, singleEnd=singleEnd)
})

setMethod("findSpliceOverlaps", c("BamFile", "ANY"),
    function(query, subject, ignore.strand=FALSE, ...,
             param=ScanBamParam(), singleEnd=TRUE)
{
    findSpliceOverlaps(.readRanges(query, param, singleEnd), subject,
                       ignore.strand, ...)
})

.readRanges <- function(bam, param, singleEnd)
{
    if (!"XS" %in% bamTag(param))
        bamTag(param) <- c(bamTag(param), "XS")
    if (singleEnd)
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
    cat("asMates:", asMates(object), "\n")
})
