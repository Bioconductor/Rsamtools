## PileupParam

.PileupParam <- setClass("PileupParam",
    representation(
        ## behavior-only policies (do not affect result schema)
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
        read_pos_breaks = "integer")) # 9

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
    function(max_depth=250, min_base_quality=0, min_mapq=13,
             min_nucleotide_depth=1, min_minor_allele_depth=0,
             distinguish_strands=TRUE, distinguish_nucleotides=TRUE,
             ignore_query_Ns=TRUE, include_deletions=TRUE,
             read_pos_breaks=integer())
{
    ## argument checking
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
    stopifnot(isTRUEorFALSE(distinguish_strands))
    stopifnot(isTRUEorFALSE(distinguish_nucleotides))
    stopifnot(isTRUEorFALSE(ignore_query_Ns))
    stopifnot(isTRUEorFALSE(include_deletions))
    if(any(is.na(read_pos_breaks)) || any(is.null(read_pos_breaks)))
        stop("'read_pos_breaks' must not contain NAs or NULLs")
    
    ## creation
    .PileupParam(max_depth=max_depth, min_base_quality=min_base_quality,
                 min_mapq=min_mapq,min_nucleotide_depth=min_nucleotide_depth,
                 min_minor_allele_depth=min_minor_allele_depth,
                 distinguish_strands=distinguish_strands,
                 distinguish_nucleotides=distinguish_nucleotides,
                 ignore_query_Ns=ignore_query_Ns,
                 include_deletions=include_deletions,
                 read_pos_breaks=read_pos_breaks)
}

.pileup <-
    function(file, index=file, ..., scanBamParam=ScanBamParam(),
             pileupParam=PileupParam())
{
    cfun <- "c_Pileup"
    if (!isOpen(file)) {
        open(file)
        on.exit(close(file))
    }
    if(bamReverseComplement(scanBamParam)) {
        warning("'reverseComplement' parameter in pileup ScanBamParam will",
                " be ignored")
        bamReverseComplement(scanBamParam) <- FALSE
    }
    schema <- .schemaBuilder(pileupParam)
    seqnamesLevels <- seqlevels(file)
    result <- .io_bam(cfun, file,
            ## space, keepFlags, isSimpleCigar extracted from 'scanBamParam'
            ## in .io_bam;
            ## remainder (...) passed to C
            bamReverseComplement(scanBamParam),
            yieldSize(file), obeyQname(file), asMates(file),
            schema, .as.list_PileupParam(pileupParam), seqnamesLevels,
            param=scanBamParam)

    which_labels <- .scanBam_extract_which_labels(scanBamParam)
    which_labels <- .make_unique(which_labels)
    run_lens <- .metacols_run_lengths(result)
    which_label <-
        rep.int(factor(which_labels, levels=which_labels), run_lens)
    
    result <- .as.data.frame_list_of_lists(result)
    result <- cbind(result, which_label) ## last col
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
                   .quickUnlist(tmp)
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

## schema
.strandHelper <- function(scanBamP) {
    flagVec <- bamFlag(scanBamP)['isMinusStrand']
    if (is.na(flagVec)) {
        "*"
    } else if (flagVec) {
        "-"
    } else {
        "+"
    }
}

.schemaBuilder <- function(pileupParam) {
    if(!inherits(pileupParam, "PileupParam")) {
        stop("'pileupParam' must inherit from 'PileupParam', got '%s'",
             class(pileupParam))
    }
    schemaDimNames <- c("strand", "nucleotide")

    strand <- ""
    if(pileupParam@distinguish_strands)
        strand <- c("+", "-")

    nucleotide <- ""
    if(pileupParam@distinguish_nucleotides) {
        nucleotide <- c("A", "C", "G", "T")
        if(!pileupParam@ignore_query_Ns)
            nucleotide <- c(nucleotide, "N")
        if(pileupParam@include_deletions)
            nucleotide <- c(nucleotide, "-")
    }

    list(schemaDimNames, list(strand, nucleotide))
}
