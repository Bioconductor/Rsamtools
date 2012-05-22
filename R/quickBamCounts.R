### =========================================================================
### quickBamCounts()
### -------------------------------------------------------------------------

.runlen <- function(x) runLength(Rle(x))

.printSummaryHeader <- function(width1, width2, width3)
{
    printHeaderLine <- function(h1, h2, h3, h4)
    {
        cat(format(h1, width=width1, justify="right"),
            " |", format(h2, width=width2, justify="right"),
            " |", format(h3, width=width3, justify="right"),
            " | ", h4, "\n", sep="")
    }
    printHeaderLine(  "group",    "nb of",  "nb of", "mean / max")
    printHeaderLine(     "of",  "records", "unique", "records per")
    printHeaderLine("records", "in group", "QNAMEs", "unique QNAME")
}

.printSummaryLine <- function(group, group_desc,
                              N_rec_per_uqname, width2, width3)
{
    N_rec <- sum(N_rec_per_uqname)
    N_uqname <- length(N_rec_per_uqname)
    if (N_uqname == 0L) {
        mean_N_rec_per_uqname <- max_N_rec_per_uqname <- NA
    } else {
        mean_N_rec_per_uqname <- mean(N_rec_per_uqname)
        max_N_rec_per_uqname <- max(N_rec_per_uqname)
    }
    cat(group_desc, " ", group,
        " |", format(N_rec, width=width2),
        " |", format(N_uqname, width=width3),
        " | ", format(mean_N_rec_per_uqname, digits=3L, width=4L),
        " / ", max_N_rec_per_uqname,
        "\n", sep="")
}

.printMainSummary <- function(N_rec_per_uqname,
                              N_1seg_rec_per_uqname,
                              N_mseg_rec_per_uqname,
                              N_first_rec_per_uqname,
                              N_last_rec_per_uqname,
                              N_other_rec_per_uqname)
{
    GROUP2DESC <- c(
        A="All records........................",
        S="  o template has single segment....",
        M="  o template has multiple segments.",
        F="      - first segment..............",
        L="      - last segment...............",
        O="      - other segment.............."
    )
    width1 <- nchar(GROUP2DESC[1L]) + 2L
    width2 <- width3 <- 9L
    .printSummaryHeader(width1, width2, width3)
    .printSummaryLine("A", GROUP2DESC[["A"]],
                      N_rec_per_uqname, width2, width3)
    .printSummaryLine("S", GROUP2DESC[["S"]],
                      N_1seg_rec_per_uqname, width2, width3)
    .printSummaryLine("M", GROUP2DESC[["M"]],
                      N_mseg_rec_per_uqname, width2, width3)
    .printSummaryLine("F", GROUP2DESC[["F"]],
                      N_first_rec_per_uqname, width2, width3)
    .printSummaryLine("L", GROUP2DESC[["L"]],
                      N_last_rec_per_uqname, width2, width3)
    .printSummaryLine("O", GROUP2DESC[["O"]],
                      N_other_rec_per_uqname, width2, width3)
    cat("\nNote that (S, M) is a partitioning of A, and (F, L, O) ",
        "is a partitioning of M.\n", "Indentation reflects this.\n", sep="")
}

.printDetailedSummary <- function(group,
                                  N_mapped_rec_per_uqname,
                                  N_primary_rec_per_uqname,
                                  N_secondary_rec_per_uqname,
                                  N_unmapped_rec_per_uqname)
{
    SUBGROUP2DESC <- c(
        `1`="  o record is mapped..............",
        `2`="      - primary alignment.........",
        `3`="      - secondary alignment.......",
        `4`="  o record is unmapped............"
    )
    width1 <- nchar(SUBGROUP2DESC[1L]) + 3L
    width2 <- width3 <- 9L
    cat("\nDetails for group ", group, ":\n", sep="")
    #.printSummaryHeader(width1, width2, width3)
    .printSummaryLine(paste0(group, "1"), SUBGROUP2DESC[["1"]],
                      N_mapped_rec_per_uqname, width2, width3)
    .printSummaryLine(paste0(group, "2"), SUBGROUP2DESC[["2"]],
                      N_primary_rec_per_uqname, width2, width3)
    .printSummaryLine(paste0(group, "3"), SUBGROUP2DESC[["3"]],
                      N_secondary_rec_per_uqname, width2, width3)
    .printSummaryLine(paste0(group, "4"), SUBGROUP2DESC[["4"]],
                      N_unmapped_rec_per_uqname, width2, width3)
}

.detailedSummary <- function(group, qnameid, flag)
{
    ## 'N_mapped_rec_per_uqname' and 'N_unmapped_rec_per_uqname':
    rec_is_unmapped <- as.logical(bamFlagAsBitMatrix(flag, "isUnmappedQuery"))
    rec_is_mapped <- !rec_is_unmapped
    N_mapped_rec_per_uqname <- .runlen(qnameid[rec_is_mapped])
    N_unmapped_rec_per_uqname <- .runlen(qnameid[rec_is_unmapped])

    ## 'N_primary_rec_per_uqname' and 'N_secondary_rec_per_uqname':
    rec_is_secondary <- as.logical(bamFlagAsBitMatrix(flag, "isNotPrimaryRead"))
    ## The notion of primary/secondary alignment only applies to mapped
    ## seqments:
    stopifnot(!any(rec_is_secondary & rec_is_unmapped))
    rec_is_primary <- rec_is_mapped & !rec_is_secondary
    N_primary_rec_per_uqname <- .runlen(qnameid[rec_is_primary])
    N_secondary_rec_per_uqname <- .runlen(qnameid[rec_is_secondary])

    .printDetailedSummary(group,
                          N_mapped_rec_per_uqname,
                          N_primary_rec_per_uqname,
                          N_secondary_rec_per_uqname,
                          N_unmapped_rec_per_uqname)
}

quickBamCounts <- function(file, index=file, param=NULL,
                           main.groups.only=FALSE)
{
    what0 <- c("qname", "flag")
    if (isSingleString(file)) {
        if (is.null(param)) {
            param <- ScanBamParam(what=what0)
        } else {
            if (length(bamWhat(param)) != 0L)
                warning("bamWhat component of supplied 'param' was ignored")
            bamWhat(param) <- what0
        }
        res <- scanBam(file, index=index, param=param)
        res0 <- res[[1L]]
        if (length(res) != 1L) {
            res0[["qname"]] <- do.call(c, lapply(res, "[[", "qname"))
            res0[["flag"]] <- do.call(c, lapply(res, "[[", "flag"))
        } 
    } else if (is.list(file) && all(what0 %in% names(file))) {
        res0 <- file
    } else {
        stop("'file' must be a single string")
    }

    if (!isTRUEorFALSE(main.groups.only))
        stop("'main.groups.only' must be TRUE or FALSE")

    ## Order records by QNAME.
    qname0 <- res0[["qname"]]
    flag0 <- res0[["flag"]]
    qnameid0 <- match(qname0, qname0)  # assign unique id to each unique QNAME
    oo <- order(qnameid0)
    qnameid <- qnameid0[oo]
    flag <- flag0[oo]

    ## 'N_rec_per_uqname':
    N_rec_per_uqname <- .runlen(qnameid)

    ## 'N_mseg_rec_per_uqname' and 'N_1seg_rec_per_uqname':
    rec_is_mseg <- as.logical(bamFlagAsBitMatrix(flag, "isPaired"))
    rec_is_1seg <- !rec_is_mseg
    N_1seg_rec_per_uqname <- .runlen(qnameid[rec_is_1seg])
    N_mseg_rec_per_uqname <- .runlen(qnameid[rec_is_mseg])
    stopifnot(identical(length(N_1seg_rec_per_uqname) +
                        length(N_mseg_rec_per_uqname),
                        length(N_rec_per_uqname)))

    ## 'N_first_rec_per_uqname' and 'N_last_rec_per_uqname' and
    ## 'N_other_rec_per_uqname':
    is_first_mate <- as.logical(bamFlagAsBitMatrix(flag, "isFirstMateRead"))
    is_second_mate <- as.logical(bamFlagAsBitMatrix(flag, "isSecondMateRead"))
    rec_is_first <- rec_is_mseg & is_first_mate & !is_second_mate
    rec_is_last <- rec_is_mseg & is_second_mate & !is_first_mate
    rec_is_other <- rec_is_mseg & (is_first_mate == is_second_mate)
    N_first_rec_per_uqname <- .runlen(qnameid[rec_is_first])
    N_last_rec_per_uqname <- .runlen(qnameid[rec_is_last])
    N_other_rec_per_uqname <- .runlen(qnameid[rec_is_other])
    stopifnot(identical(sum(N_first_rec_per_uqname) +
                        sum(N_last_rec_per_uqname) +
                        sum(N_other_rec_per_uqname),
                        sum(N_mseg_rec_per_uqname)))

    .printMainSummary(N_rec_per_uqname,
                      N_1seg_rec_per_uqname,
                      N_mseg_rec_per_uqname,
                      N_first_rec_per_uqname,
                      N_last_rec_per_uqname,
                      N_other_rec_per_uqname)

    if (main.groups.only)
        return(invisible(NULL))
    if (length(N_1seg_rec_per_uqname) != 0L &&
        length(N_mseg_rec_per_uqname != 0L))
        .detailedSummary("A", qnameid, flag)
    if (length(N_1seg_rec_per_uqname) != 0L)
        .detailedSummary("S", qnameid[rec_is_1seg], flag[rec_is_1seg])
    if (length(N_mseg_rec_per_uqname) != 0L)
        .detailedSummary("M", qnameid[rec_is_mseg], flag[rec_is_mseg])
    if (length(N_first_rec_per_uqname) != 0L)
        .detailedSummary("F", qnameid[rec_is_first], flag[rec_is_first])
    if (length(N_last_rec_per_uqname) != 0L)
        .detailedSummary("L", qnameid[rec_is_last], flag[rec_is_last])
    if (length(N_other_rec_per_uqname) != 0L)
        .detailedSummary("O", qnameid[rec_is_other], flag[rec_is_other])
    return(invisible(NULL))
}

