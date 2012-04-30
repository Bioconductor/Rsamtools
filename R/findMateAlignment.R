### =========================================================================
### findMateAlignment()
### -------------------------------------------------------------------------
### For each element in GappedAlignments object 'x', finds its mate in
### GappedAlignments object 'y'.
### 'x[i1]' and 'y[i2]' are considered mates iff:
###   (A) names(x[i1]) == names(y[i2])
###   (B) elementMetadata(x[i1])$mrnm == seqnames(y[i2]) &
###       elementMetadata(y[i2])$mrnm == seqnames(x[i1])
###   (C) elementMetadata(x[i1])$mpos == start(y[i2]) &
###       elementMetadata(y[i2])$mpos == start(x[i1])
###   (D) isMateMinusStrand(x[i1]) == isMinusStrand(y[i2]) &
###       isMateMinusStrand(y[i2]) == isMinusStrand(x[i1])
###   (E) isFirstSegment(x[i1]) & isLastSegment(y[i2]) |
###       isFirstSegment(y[i2]) & isLastSegment(x[i1])

.checkElementMetadata <- function(arg, argname)
{
    if (!is(arg, "GappedAlignments"))
        stop("'", argname, "' must be a GappedAlignments object")
    if (is.null(names(arg)))
        stop("'", argname, "' must have names")
    eltmetadata <- elementMetadata(arg)
    REQUIRED_COLNAMES <- c("flag", "mrnm", "mpos")
    if (!all(REQUIRED_COLNAMES %in% colnames(eltmetadata))) {
        colnames_in1string <- paste("\"", REQUIRED_COLNAMES, "\"", sep="",
                                    collapse=", ")
        stop("required columns in 'elementMetadata(", argname, ")': ",
             colnames_in1string)
    }
    eltmetadata
}

### 'names', 'flagbits', 'mrnm', and 'mpos', must all come from the same
### GappedAlignments object x.
### 'names': names(x).
### 'flagbits': integer matrix (of 0's and 1's) obtained with
###     bamFlagAsBitMatrix(elementMetadata(x)$flag,
###                        bitnames=.MATING_FLAG_BITNAMES)
### 'mrnm': character vector or factor obtained with elementMetadata(x)$mrnm
### 'mpos': integer vector obtained with elementMetadata(x)$mpos
### Returns 'names' with NAs injected at positions corresponding to alignments
### that satisfy at least one of following conditions:
###     1. Bit 0x1 (isPaired) is 0
###     2. Read is neither first or last mate
###     3. Bit 0x8 (hasUnmappedMate) is 1
###     4. 'mrnm' is NA (i.e. RNEXT = '*')
###     5. 'mpos' is NA (i.e. PNEXT = 0)
### My understanding of the SAM Spec is that 3., 4. and 5. should happen
### simultaneously even though the Spec don't clearly state this.

.MATING_FLAG_BITNAMES <- c("isPaired", "hasUnmappedMate",
                           "isFirstMateRead", "isSecondMateRead")

.makeGappedAlignmentsGNames <- function(names, flagbits, mrnm, mpos)
{
    is_paired <- flagbits[ , "isPaired"]
    is_first <- flagbits[ , "isFirstMateRead"]
    is_last <- flagbits[ , "isSecondMateRead"]
    has_unmappedmate <- flagbits[ , "hasUnmappedMate"]
    alter_idx <- which(!is_paired |
                       is_first == is_last |
                       has_unmappedmate |
                       is.na(mrnm) |
                       is.na(mpos))
    names[alter_idx] <- NA_integer_
    names
}

### Puts NAs last.
.getCharacterOrderAndGroupSizes <- function(x)
{
    x2 <- match(x, x,
                nomatch=.Machine$integer.max,
                incomparables=NA_character_)
    xo <- IRanges:::orderInteger(x2)
    ox2 <- Rle(x2[xo])
    group.sizes <- runLength(ox2)
    ngroup <- length(group.sizes)
    if (ngroup != 0L && runValue(ox2)[ngroup] == .Machine$integer.max)
        group.sizes <- group.sizes[-ngroup]
    list(xo=xo, group.sizes=group.sizes)
}

### Assumes that GappedAlignments object 'x':
###   (1) has clean group names i.e. .makeGappedAlignmentsGNames() would not
###       inject any NA in its names,
###   (2) is already ordered by group names (which are the same as its names).
### 'group.sizes' must be the same as 'runLength(Rle(names(x)))'.
.findMateAlignmentInChunk <- function(group.sizes,
                                      x_mrnm, x_seqnames,
                                      x_mpos, x_start,
                                      x_is_mate_minus, x_is_minus,
                                      x_is_first)
{
    hits <- IRanges:::makeAllGroupInnerHits(group.sizes, hit.type=1L)
    x_hits <- queryHits(hits)
    y_hits <- subjectHits(hits)

    ## Keep hits that satisfy conditions (B), (C), (D) and (E).
    valid_hits <- x_mrnm[x_hits] == x_seqnames[y_hits] &  # (B)
                  x_mrnm[y_hits] == x_seqnames[x_hits] &  # (B)
                  x_mpos[x_hits] == x_start[y_hits] &  # (C)
                  x_mpos[y_hits] == x_start[x_hits] &  # (C)
                  x_is_mate_minus[x_hits] == x_is_minus[y_hits] &  # (D)
                  x_is_mate_minus[y_hits] == x_is_minus[x_hits] &  # (D)
                  x_is_first[x_hits] != x_is_first[y_hits]  # (E)
    x_hits <- x_hits[valid_hits]
    y_hits <- y_hits[valid_hits]

    tmp <- x_hits
    x_hits <- c(x_hits, y_hits)
    y_hits <- c(y_hits, tmp)

    ans <- rep.int(NA_integer_, length(x_start))
    ans[x_hits] <- y_hits
    is_dup <- duplicated(x_hits)
    ans[x_hits[is_dup]] <- -1L
    ans
}

### Takes about 6 sec and 274MB of RAM to mate 1 million alignments,
### and about 26.3 sec and 1022MB of RAM to mate 5 million alignments.
findMateAlignment <- function(x, verbose=FALSE)
{
    if (!isTRUEorFALSE(verbose))
        stop("'verbose' must be TRUE or FALSE")
    x_names <- names(x)
    if (is.null(x_names))
        stop("'x' must have names")
    x_eltmetadata <- .checkElementMetadata(x, "x")
    x_flag <- x_eltmetadata$flag
    bitnames <- c(.MATING_FLAG_BITNAMES, "isMinusStrand", "isMateMinusStrand")
    x_flagbits <- bamFlagAsBitMatrix(x_flag, bitnames=bitnames)
    x_mrnm <- as.character(x_eltmetadata$mrnm)
    x_mpos <- x_eltmetadata$mpos
    x_gnames <- .makeGappedAlignmentsGNames(x_names, x_flagbits,
                                            x_mrnm, x_mpos)
    x_seqnames <- as.character(seqnames(x))
    x_start <- start(x)
    x_is_mate_minus <- x_flagbits[ , "isMateMinusStrand"]
    x_is_minus <- x_flagbits[ , "isMinusStrand"]
    x_is_first <- x_flagbits[ , "isFirstMateRead"]

    xo_and_GS <- .getCharacterOrderAndGroupSizes(x_gnames)
    xo <- xo_and_GS$xo
    GS <- xo_and_GS$group.sizes
    ans <- rep.int(NA_integer_, length(x_gnames))
    NGROUP_BY_CHUNK <- 25000L
    chunk.GIDX <- seq_len(NGROUP_BY_CHUNK)
    chunk.offset <- 0L
    while (TRUE) {
        chunk.GIDX <- chunk.GIDX[chunk.GIDX <= length(GS)]
        if (length(chunk.GIDX) == 0L)
            break
        chunk.GS <- GS[chunk.GIDX]
        chunk.length <- sum(chunk.GS)
        chunk.idx <- xo[chunk.offset + seq_len(chunk.length)]
        chunk.x_mrnm <- x_mrnm[chunk.idx]
        chunk.x_seqnames <- x_seqnames[chunk.idx]
        chunk.x_mpos <- x_mpos[chunk.idx]
        chunk.x_start <- x_start[chunk.idx]
        chunk.x_is_mate_minus <- x_is_mate_minus[chunk.idx]
        chunk.x_is_minus <- x_is_minus[chunk.idx]
        chunk.x_is_first <- x_is_first[chunk.idx]
        if (verbose)
            message("Finding mates in chunk of ", chunk.length,
                    " alignments ... ", appendLF=FALSE)
        chunk.ans <- .findMateAlignmentInChunk(chunk.GS,
                                               chunk.x_mrnm,
                                               chunk.x_seqnames,
                                               chunk.x_mpos,
                                               chunk.x_start,
                                               chunk.x_is_mate_minus,
                                               chunk.x_is_minus,
                                               chunk.x_is_first)
        have_more_than_1_mate <- which(chunk.ans < 0L)
        if (length(have_more_than_1_mate) != 0L) {
            morethan1mate_idx <- chunk.idx[have_more_than_1_mate]
            cat("\n!! findMateAlignment() found more than 1 mate in 'x' ",
                "for elements: ",
                paste(morethan1mate_idx, collapse=", "),
                ".\n!! Details:\n!! ", sep="")
            GenomicRanges:::showGappedAlignments(x[morethan1mate_idx],
                                                 margin="!! ",
                                                 with.classinfo=TRUE,
                                                 print.seqlengths=FALSE)
            cat("!! ==> won't assign a mate to those elements!\n")
            chunk.ans[have_more_than_1_mate] <- NA_integer_
        }
        if (verbose)
            message("OK")
        ans[chunk.idx] <- chunk.idx[chunk.ans]
        chunk.GIDX <- chunk.GIDX + NGROUP_BY_CHUNK
        chunk.offset <- chunk.offset + chunk.length
    }
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### findMateAlignment2().
###

### .findMatches() is the same as match() except that it returns *all*
### the matches (in a Hits object, ordered by queryHits first, then by
### subjectHits).
### TODO: Make findMatches() an S4 generic function with at least a method for
### vectors. Like findOverlaps(), findMatches() could support the 'select' arg
### (but with supported values "all", "first" and "last" only, no need for
### "arbitrary") so that when used with 'select="first"', it would be
### equivalent to match(). This stuff would go in IRanges.
.findMatches <- function(query, subject, incomparables=NULL)
{
    if (!is.vector(query) || !is.vector(subject))
        stop("'query' and 'subject' must be vectors")
    if (class(query) != class(subject))
        stop("'query' and 'subject' must be vectors of the same class")
    if (!is.null(incomparables) && !(is.vector(incomparables) &&
                                     class(incomparables) == class(query)))
        stop("'incomparables' must be NULL or a vector ",
             "of the same class as 'query' and 'subject'")
    m0 <- match(query, subject, incomparables=incomparables)
    query_hits0 <- which(!is.na(m0))
    if (length(query_hits0) == 0L) {
        query_hits <- subject_hits <- integer(0)
    } else {
        subject_hits0 <- m0[query_hits0]
        subject_low2high <- IRanges:::.makeLow2highFromHigh2low(
                                high2low(subject))
        extra_hits <- subject_low2high[subject_hits0]
        query_nhits <- 1L + elementLengths(extra_hits)
        query_hits <- rep.int(query_hits0, query_nhits)
        subject_hits <- integer(length(query_hits))
        idx0 <- cumsum(c(1L, query_nhits[-length(query_nhits)]))
        subject_hits[idx0] <- m0[query_hits0]
        subject_hits[-idx0] <- unlist(extra_hits,
                                      recursive=FALSE, use.names=FALSE)
    }
    new2("Hits", queryHits=query_hits, subjectHits=subject_hits,
                 queryLength=length(query), subjectLength=length(subject),
                 check=FALSE)
}

### Use to find self matches in 'x'. Twice faster than
### 'findMatches(x, x, incomparables=NA_character_)' and uses
### twice less memory.
.findSelfMatches.character <- function(x)
{
    xo_and_GS <- .getCharacterOrderAndGroupSizes(x)
    xo <- xo_and_GS$xo
    GS <- xo_and_GS$group.sizes
    ans <- IRanges:::makeAllGroupInnerHits(GS, hit.type=1L)
    ans@queryHits <- xo[ans@queryHits]
    ans@subjectHits <- xo[ans@subjectHits]
    ans@queryLength <- ans@subjectLength <- length(x)
    ans
}

### Takes about 5.47 sec and 295MB of RAM to mate 1 million alignments,
### and about 32.26 sec and 1330MB of RAM to mate 5 million alignments.
findMateAlignment2 <- function(x, y=NULL)
{
    x_names <- names(x)
    if (is.null(x_names))
        stop("'x' must have names")
    x_eltmetadata <- .checkElementMetadata(x, "x")
    x_flag <- x_eltmetadata$flag
    bitnames <- c(.MATING_FLAG_BITNAMES, "isMinusStrand", "isMateMinusStrand")
    x_flagbits <- bamFlagAsBitMatrix(x_flag, bitnames=bitnames)
    x_mrnm <- as.character(x_eltmetadata$mrnm)
    x_mpos <- x_eltmetadata$mpos
    x_gnames <- .makeGappedAlignmentsGNames(x_names, x_flagbits,
                                            x_mrnm, x_mpos)
    x_seqnames <- as.character(seqnames(x))
    x_start <- start(x)
    x_is_mate_minus <- x_flagbits[ , "isMateMinusStrand"]
    x_is_minus <- x_flagbits[ , "isMinusStrand"]
    x_is_first <- x_flagbits[ , "isFirstMateRead"]

    if (is.null(y)) {
        y_mrnm <- x_mrnm
        y_seqnames <- x_seqnames
        y_mpos <- x_mpos
        y_start <- x_start
        y_is_mate_minus <- x_is_mate_minus
        y_is_minus <- x_is_minus
        y_is_first <- x_is_first

        hits <- .findSelfMatches.character(x_gnames)
    } else {
        y_names <- names(y)
        if (is.null(y_names))
            stop("'y' must have names")
        y_eltmetadata <- .checkElementMetadata(y, "y")
        y_flag <- y_eltmetadata$flag
        y_flagbits <- bamFlagAsBitMatrix(y_flag, bitnames=bitnames)
        y_mrnm <- as.character(y_eltmetadata$mrnm)
        y_mpos <- y_eltmetadata$mpos
        y_gnames <- .makeGappedAlignmentsGNames(y_names, y_flagbits,
                                                y_mrnm, y_mpos)
        y_seqnames <- as.character(seqnames(y))
        y_start <- start(y)
        y_is_mate_minus <- y_flagbits[ , "isMateMinusStrand"]
        y_is_minus <- y_flagbits[ , "isMinusStrand"]
        y_is_first <- y_flagbits[ , "isFirstMateRead"]

        hits <- .findMatches(x_gnames, y_gnames, incomparables=NA_character_)
    }

    x_hits <- queryHits(hits)
    y_hits <- subjectHits(hits)

    ## Keep hits that satisfy conditions (B), (C), (D) and (E).
    valid_hits <- x_mrnm[x_hits] == y_seqnames[y_hits] &  # (B)
                  y_mrnm[y_hits] == x_seqnames[x_hits] &  # (B)
                  x_mpos[x_hits] == y_start[y_hits] &  # (C)
                  y_mpos[y_hits] == x_start[x_hits] &  # (C)
                  x_is_mate_minus[x_hits] == y_is_minus[y_hits] &  # (D)
                  y_is_mate_minus[y_hits] == x_is_minus[x_hits] &  # (D)
                  x_is_first[x_hits] != y_is_first[y_hits]  # (E)
    x_hits <- x_hits[valid_hits]
    y_hits <- y_hits[valid_hits]

    if (is.null(y)) {
        tmp <- x_hits
        x_hits <- c(x_hits, y_hits)
        y_hits <- c(y_hits, tmp)
    }
    is_dup <- duplicated(x_hits)
    if (any(is_dup)) {
        have_more_than_one_mate <- unique(x_hits[is_dup])
        stop("more than 1 mate found for elements ",
             paste(have_more_than_one_mate, collapse=", "))
    }
    ans <- rep.int(NA_integer_, length(x))
    ans[x_hits] <- y_hits
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeGappedAlignmentPairs().
###

### TODO: Make isFirstSegment() an S4 generic function with methods for
### matrices, integer vectors, and GappedAlignments objects. Put this with the
### flag utils in Rsamtools.
.isFirstSegment.matrix <- function(x)
{
    is_paired <- as.logical(x[ , "isPaired"])
    is_first0 <- as.logical(x[ , "isFirstMateRead"])
    is_last0 <- as.logical(x[ , "isSecondMateRead"])
    ## According to SAM Spec, bits 0x40 (isFirstMateRead) and 0x80
    ## (isSecondMateRead) can both be set or unset, even when bit 0x1
    ## (isPaired) is set. However we are not interested in those situations
    ## (which have a special meaning).
    is_paired & is_first0 & (!is_last0)
}

.isFirstSegment.integer <- function(flag)
{
    bitnames <- c("isPaired", "isFirstMateRead", "isSecondMateRead")
    .isFirstSegment.matrix(bamFlagAsBitMatrix(flag, bitnames=bitnames))
}

.isFirstSegment.GappedAlignments <- function(x)
    .isFirstSegment.integer(elementMetadata(x)$flag)

### TODO: Make isLastSegment() an S4 generic function with methods for
### matrices, integer vectors, and GappedAlignments objects. Put this with the
### flag utils in Rsamtools.
.isLastSegment.matrix <- function(x)
{
    is_paired <- as.logical(x[ , "isPaired"])
    is_first0 <- as.logical(x[ , "isFirstMateRead"])
    is_last0 <- as.logical(x[ , "isSecondMateRead"])
    ## According to SAM Spec, bits 0x40 (isFirstMateRead) and 0x80
    ## (isSecondMateRead) can both be set or unset, even when bit 0x1
    ## (isPaired) is set. However we are not interested in those situations
    ## (which have a special meaning).
    is_paired & is_last0 & (!is_first0)
}

.isLastSegment.integer <- function(flag)
{
    bitnames <- c("isPaired", "isFirstMateRead", "isSecondMateRead")
    .isLastSegment.matrix(bamFlagAsBitMatrix(flag, bitnames=bitnames))
}

.isLastSegment.GappedAlignments <- function(x)
    .isLastSegment.integer(elementMetadata(x)$flag)

### 'x' must be a GappedAlignments objects.
makeGappedAlignmentPairs <- function(x, use.names=FALSE, keep.cols=NULL)
{
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    if (!is.null(keep.cols)) {
        if (!is.character(keep.cols) || any(is.na(keep.cols)))
            stop("'keep.cols' must be a character vector with no NAs")
        if (!all(keep.cols %in% colnames(elementMetadata(x))))
            stop("'keep.cols' must be a subset ",
                 "of 'colnames(elementMetadata(x))'")
    }
    mate <- findMateAlignment(x)
    x_is_first <- .isFirstSegment.GappedAlignments(x)
    x_is_last <- .isLastSegment.GappedAlignments(x)
    first_idx <- which(!is.na(mate) & x_is_first)
    last_idx <- mate[first_idx]

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

    ## Check the 0x2 bit (isProperPair).
    x_is_proper <- as.logical(bamFlagAsBitMatrix(elementMetadata(x)$flag,
                                                 bitnames="isProperPair"))
    first_is_proper <- x_is_proper[first_idx]
    last_is_proper <- x_is_proper[last_idx]
    if (!identical(first_is_proper, last_is_proper))
        stop("for some pairs, the 2 mates don't have ",
             "the same \"isProperPair\" bit.\n",
             "  Maybe the BAM file they are coming from was corrupted or ",
             "generated by\n  an aligner that doesn't follow the SAM Spec?")

    ## Drop pairs with discordant seqnames or strand.
    idx_is_discordant <- (as.character(seqnames(x)[first_idx]) !=
                          as.character(seqnames(x)[last_idx])) |
                         (as.character(strand(x)[first_idx]) ==
                          as.character(strand(x)[last_idx]))
    if (any(idx_is_discordant) != 0L) {
        nb_discordant_proper <- sum(first_is_proper[idx_is_discordant])
        if (nb_discordant_proper != 0L) {
            ratio <- 100.0 * nb_discordant_proper / sum(idx_is_discordant)
            warning(ratio, "% of the pairs with discordant seqnames or ",
                    "strand were flagged\n",
                    "  as proper pairs by the aligner. Dropping them anyway.")
        }
        keep <- -which(idx_is_discordant)
        first_idx <- first_idx[keep]
        last_idx <- last_idx[keep]
        first_is_proper <- first_is_proper[keep]
    }

    ## The big split!
    first <- x[first_idx]
    last <- x[last_idx]
    ans_names <- NULL
    if (use.names)
        ans_names <- names(first)
    names(first) <- names(last) <- NULL
    elementMetadata(first) <- elementMetadata(first)[keep.cols]
    elementMetadata(last) <- elementMetadata(last)[keep.cols]
    GappedAlignmentPairs(first, last, first_is_proper, names=ans_names)
}

