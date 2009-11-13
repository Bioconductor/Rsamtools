###
cigarToReadWidth <- function(cigar, after.hard.clipping=FALSE)
{
    if (!is.character(cigar)) {
        if (!is.factor(cigar) || !is.character(levels(cigar)))
            stop("'cigar' must be a character vector/factor")
        cigar <- as.vector(cigar)
    }
    if (!isTRUEorFALSE(after.hard.clipping))
        stop("'after.hard.clipping' must be TRUE or FALSE")
    .Call(".cigar_to_read_width",
          cigar, after.hard.clipping,
          PACKAGE="Rsamtools")
}

cigarToIRanges <- function(cigar, drop.D.ranges=FALSE, merge.ranges=TRUE)
{
    if (is.factor(cigar) && is.character(levels(cigar)))
        cigar <- as.vector(cigar)
    if (!isSingleString(cigar))
        stop("'cigar' must be a single string")
    if (!isTRUEorFALSE(drop.D.ranges))
        stop("'drop.D.ranges' must be TRUE or FALSE")
    if (!isTRUEorFALSE(merge.ranges))
        stop("'merge.ranges' must be TRUE or FALSE")
    .Call(".cigar_to_IRanges", cigar, drop.D.ranges, merge.ranges,
          PACKAGE="Rsamtools")
}

cigarToIRangesList <- function(cigar, rname, pos, flag=NULL,
                               drop.D.ranges=FALSE, merge.ranges=TRUE)
{
    if (!is.character(cigar)) {
        if (!is.factor(cigar) || !is.character(levels(cigar)))
            stop("'cigar' must be a character vector/factor")
        cigar <- as.vector(cigar)
    }
    if (!is.factor(rname) || !is.character(levels(rname))) {
        if (!is.character(rname))
            stop("'rname' must be a character vector/factor")
        rname <- as.factor(rname)
    }
    if (!is.numeric(pos))
        stop("'pos' must be a vector of integers")
    if (!is.integer(pos))
        pos <- as.integer(pos)
    if (length(cigar) != length(rname) || length(cigar) != length(pos))
        stop("'cigar', 'rname' and 'pos' must have the same length")
    if (!is.null(flag)) {
        if (!is.numeric(flag))
            stop("'flag' must be NULL or a vector of integers")
        if (!is.integer(flag))
            flag <- as.integer(flag)
        if (length(cigar) != length(flag))
            stop("'cigar' and 'flag' must have the same length")
    }
    if (!isTRUEorFALSE(drop.D.ranges))
        stop("'drop.D.ranges' must be TRUE or FALSE")
    if (!isTRUEorFALSE(merge.ranges))
        stop("'merge.ranges' must be TRUE or FALSE")
    C_ans <- .Call(".cigar_to_list_of_IRanges",
                   cigar, rname, pos, flag, drop.D.ranges, merge.ranges,
                   PACKAGE="Rsamtools")
    if (length(C_ans) < 200L)
        IRangesList(C_ans, compress=FALSE)
    else
        IRangesList(C_ans, compress=TRUE)
}

cigarToCigarTable <- function(cigar) {
    if (is.character(cigar))
        cigar <- factor(cigar)
    else if (!is.factor(cigar))
        stop("'cigar' must be a character vector/factor")
    basicTable <- table(cigar)
    tableOrder <- order(basicTable, decreasing=TRUE)
    cigar <- factor(cigar, levels = levels(cigar)[tableOrder])
    basicTable <- basicTable[tableOrder]
    cigarValues <-
      CharacterList(lapply(strsplit(levels(cigar), "[0-9]+"), "[", -1))
    cigarLengths <- IntegerList(strsplit(levels(cigar), "[A-Za-z]+"))
    DataFrame(cigar =
              IRanges:::newCompressedList("CompressedRleList",
                                   Rle(unlist(cigarValues, use.names=FALSE),
                                       unlist(cigarLengths, use.names=FALSE)),
                                   cumsum(unlist(lapply(cigarLengths, sum)))),
              count = as.integer(basicTable))
}

summarizeCigarTable <- function(x) {
    alignedCharacters <-
      table(rep.int(elementLengths(x[["cigar"]]), x[["count"]]),
            rep.int(viewSums(Views(unlist(x[["cigar"]] == "M"),
                                   as(x[["cigar"]]@partitioning, "IRanges"))) ==
                                   elementLengths(x[["cigar"]]),
                    x[["count"]]))
    tabledAlignedCharacters <- as(rev(colSums(alignedCharacters)), "integer")
    names(tabledAlignedCharacters) <-
      unname(c("TRUE" = "AllAligned",
               "FALSE" = "SomeNonAligned")[names(tabledAlignedCharacters)])

    indelHits <-
      rbind(data.frame(subject =
                       subjectHits(findOverlaps(IRanges(unlist(x[["cigar"]] == "D")),
                                                x[["cigar"]]@partitioning)),
                       type = factor("D", levels = c("D", "I"))),
            data.frame(subject =
                       subjectHits(findOverlaps(IRanges(unlist(x[["cigar"]] == "I")),
                                                x[["cigar"]]@partitioning)),
                       type = factor("I", levels = c("D", "I"))))
    tabledIndelHits <- table(indelHits[,1], indelHits[,2])
    tabledIndelHits <-
      tabledIndelHits[rep.int(seq_len(nrow(tabledIndelHits)),
                              x[["count"]][as.integer(rownames(tabledIndelHits))]),
                      , drop = FALSE]
    tabledIndelHits <-
      as(table(tabledIndelHits[,"D"], tabledIndelHits[,"I"]), "matrix")
    rownames(tabledIndelHits) <- paste("D", rownames(tabledIndelHits), sep = "")
    colnames(tabledIndelHits) <- paste("I", colnames(tabledIndelHits), sep = "")
    if (!("D0" %in% rownames(tabledIndelHits)))
        tabledIndelHits <-
          rbind("D0" = rep.int(0L, ncol(tabledIndelHits)), tabledIndelHits)
    if (!("I0" %in% colnames(tabledIndelHits)))
        tabledIndelHits <-
          cbind("I0" = rep.int(0L, nrow(tabledIndelHits)), tabledIndelHits)
    tabledIndelHits["D0", "I0"] <- nrow(x) - sum(tabledIndelHits[-1])

    list("AlignedCharacters" = tabledAlignedCharacters,
         "Indels" = tabledIndelHits)
}
