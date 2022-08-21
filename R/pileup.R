## PileupParam

.PileupParam <- setClass("PileupParam",
    representation(
        ## behavior-only policies
        max_depth = "integer", #pileupParams elt 0 (in C code)
        min_base_quality = "integer", # 1
        min_mapq = "integer", # 2
        min_nucleotide_depth = "integer", # 3
        min_minor_allele_depth = "integer", # 4
        
        ## behavior+structure policies
        distinguish_strands = "logical", # 5
        distinguish_nucleotides = "logical", # 6
        ignore_query_Ns = "logical", # 7
        include_deletions="logical", # 8
        include_insertions="logical", # 9
        left_bins = "numeric", # 10
        query_bins="numeric")) # 11

setMethod(show, "PileupParam", function(object) {
    cat("class: ", class(object), "\n")
    values <- sapply(slotNames(object), slot, object=object)
    info <- paste(slotNames(object), values, sep=": ", collapse="; ")
    cat(strwrap(info, exdent=2), sep="\n")
})

.as.list_PileupParam <- function(x, ...) {
    slotnames <- slotNames(x)
    names(slotnames) <- slotnames
    lapply(slotnames, slot, object=x)
}

PileupParam <-
    function(max_depth=250, min_base_quality=13, min_mapq=0,
             min_nucleotide_depth=1, min_minor_allele_depth=0,
             distinguish_strands=TRUE, distinguish_nucleotides=TRUE,
             ignore_query_Ns=TRUE, include_deletions=TRUE,
             include_insertions=FALSE, left_bins=NULL,
             query_bins=NULL, cycle_bins=NULL)
{
    ## argument checking
    if(!is.null(cycle_bins)) {
        .Deprecated("cycle_bins", "Rsamtools",
                    paste("'cycle_bins' deprecated; rename 'cycle_bins'",
                          "to 'left_bins' (identical behvaior)"))
        left_bins <- cycle_bins
    }
    if(!is.null(left_bins) && !is.null(query_bins))
        stop("only one of 'left_bins', 'query_bins', and 'cycle_bins' allowed")
    if(is.null(left_bins)) left_bins <- numeric()
    if(is.null(query_bins)) query_bins <- numeric()
    ## invariant:
    ## - left_bins & query_bins length 0
    ## - one of left_bins or query_bins length > 0, but the other length 0

    stopifnot(isSingleNumber(max_depth))
    stopifnot(isSingleNumber(min_base_quality))
    stopifnot(isSingleNumber(min_mapq))
    stopifnot(isSingleNumber(min_nucleotide_depth))
    stopifnot(isSingleNumber(min_minor_allele_depth))
    max_depth <- as.integer(max_depth)
    min_base_quality <- as.integer(min_base_quality)
    min_mapq <- as.integer(min_mapq)
    min_nucleotide_depth <- as.integer(min_nucleotide_depth)
    min_minor_allele_depth <- as.integer(min_minor_allele_depth)
    left_bins <- .preprocess_bins(left_bins)
    query_bins <- .preprocess_bins(query_bins)

    stopifnot(isTRUEorFALSE(distinguish_strands))
    stopifnot(isTRUEorFALSE(distinguish_nucleotides))
    stopifnot(isTRUEorFALSE(ignore_query_Ns))
    stopifnot(isTRUEorFALSE(include_deletions))
    stopifnot(isTRUEorFALSE(include_insertions))
    
    ## creation
    .PileupParam(max_depth=max_depth, min_base_quality=min_base_quality,
                 min_mapq=min_mapq,min_nucleotide_depth=min_nucleotide_depth,
                 min_minor_allele_depth=min_minor_allele_depth,
                 distinguish_strands=distinguish_strands,
                 distinguish_nucleotides=distinguish_nucleotides,
                 ignore_query_Ns=ignore_query_Ns,
                 include_deletions=include_deletions,
                 include_insertions=include_insertions, left_bins=left_bins,
                 query_bins=query_bins)
}

.pileup <-
    function(file, index=file, ..., scanBamParam=ScanBamParam(),
             pileupParam=PileupParam())
{
    if (!isOpen(file)) {
        open(file)
        on.exit(close(file))
    }
    if(bamReverseComplement(scanBamParam)) {
        warning("'reverseComplement' parameter in pileup ScanBamParam will",
                " be ignored")
        bamReverseComplement(scanBamParam) <- FALSE
    }
    result <- .io_bam(.c_Pileup, file,
            ## space, keepFlags, isSimpleCigar extracted from 'scanBamParam'
            ## in .io_bam;
            ## remainder (...) passed to C
            bamReverseComplement(scanBamParam),
            yieldSize(file), obeyQname(file), asMates(file),
            qnamePrefixEnd(file), qnameSuffixStart(file),
            .as.list_PileupParam(pileupParam), param=scanBamParam)
    ##browser()

    which_labels <- .scanBam_extract_which_labels(scanBamParam)
    which_labels <- .make_unique(which_labels)
    run_lens <- .metacols_run_lengths(result)
    which_label <-
        rep.int(factor(which_labels, levels=which_labels), run_lens)

    ## no-op if no bins
    if(length(pileupParam@left_bins) > 0)
        result <- .apply_bin_levels(result, pileupParam@left_bins)
    else if(length(pileupParam@query_bins) > 0)
        result <- .apply_bin_levels(result, pileupParam@query_bins)

    result <- .as.data.frame_list_of_lists(result)

    ## wait to rename column until after converted to data.frame
    if(length(pileupParam@left_bins) > 0) {
        if(! "bin" %in% names(result))
            stop("internal: expected a 'bin' column to rename 'left_bin'")
        names(result)[names(result) == "bin"] <- "left_bin"
    } else if(length(pileupParam@query_bins) > 0) {
        if(! "bin" %in% names(result))
            stop("internal: expected a 'bin' column to rename 'query_bin'")
        names(result)[names(result) == "bin"] <- "query_bin"
    }

    if(length(bamWhich(scanBamParam)) != 0L) ## if no space arg
        result <- cbind(result, which_label) ## last col
    result
}

.apply_bin_levels <- function(result, bins) {
    ## no-op if user didn't ask for bins
    if(length(bins) > 0L) {
        bin_levels <- .bin_levels(bins)
        for(i in seq_along(result))
            result[[i]]$bin <- .as.factor_bin(result[[i]]$bin, bin_levels)
    }
    result
}

.as.factor_bin <- function(bin, bin_levels) {
    structure(bin, levels=bin_levels, class="factor")
    ## equivalent:
    ## attributes(cycle_bin) <- list(levels=cycle_bin_levels, class="factor")
    ## cycle_bin
}

.bin_levels <- function(bins) {
    bins[bins == .Machine$integer.max] <- Inf
    bins[bins == -.Machine$integer.max] <- -Inf
    levels(cut(0, bins))
}

## return value: numeric vector of increasing values that contains
## only integers and +/-Inf
.preprocess_bins <- function(bins) {
    if(length(bins) != 0L) {
        if(any(is.na(bins)) || any(is.null(bins)) || any(is.nan(bins)))
            stop("bin args must not contain NAs, NULLs, or NaNs")
        if(length(bins) == 1L)
            stop("bins args must have 0 or >1 elements")
        if(any(bins < 0L) && any(bins > 0L))
            stop("bin args values must all have the same sign (or be Inf)")
        if(any(bins == 0L) && any(bins < 0L))
            stop("'0' not allowed when specifying reverse bins; try '-1'?")
        ## invariant: only contains integers and +/-Inf
        bins[!is.finite(bins) & bins > 0] <- .Machine$integer.max ## Inf to max_int
        bins[!is.finite(bins) & bins < 0] <- -.Machine$integer.max ## -Inf to -max_int
        bins <- as.integer(bins)
        if(any(duplicated(bins)))
           stop("bin args must not contain duplicate values")
        bins <- sort(bins)
    }
    bins
}

.make_unique <- function(group_names) {
    if(anyDuplicated(group_names))
        group_names <- paste0(group_names, ".", seq_along(group_names))
    else
        group_names
}

.metacols_run_lengths <- function(x) {
    if(length(x) == 0L)
        run_lens <- integer(0)
    else
        run_lens <- sapply(x, function(x_elt) length(x_elt[[1L]]))
}

.as.data.frame_list_of_lists <-function(x) {
    if(length(x) < 1L) {
        stop("'x' must have length > 0, got '%d'", length(x))
    }
    ans_colnames <- names(x[[1L]])
    ans <-
        lapply(setNames(ans_colnames, ans_colnames),
               function(nm) {
                   tmp <- lapply(x, "[[", nm)
                   S4Vectors:::quick_unlist(tmp)
               })
    ans <- data.frame(ans, stringsAsFactors=FALSE)
}

setMethod("pileup", "character",
    function(file, index=file, ..., scanBamParam=ScanBamParam(),
             pileupParam=PileupParam())
{
    stopifnot(length(file) == 1L)
    bf <- BamFile(file, index=index)
    .pileup(bf, scanBamParam=scanBamParam, pileupParam=pileupParam)
})

setMethod("pileup", "BamFile", .pileup)

.pileupWhat <- function(pileupParam) {
    result <- c("pos",
      if (pileupParam@distinguish_strands) "strand" else NULL,
      if (pileupParam@distinguish_nucleotides) "nucleotide" else NULL,
      "count")
    setNames(result, result)
}
