### =========================================================================
### findMateAlignment()
### -------------------------------------------------------------------------
###
### For each element in GappedAlignments object 'x', finds its mate in
### GappedAlignments object 'y'.
###
### Alignments 'x[i1]' and 'y[i2]' are considered mates iff:
###
###   (A) names(x[i1]) == names(y[i2])
###
###   (B) mcols(x[i1])$mrnm == seqnames(y[i2]) &
###       mcols(y[i2])$mrnm == seqnames(x[i1])
###
###   (C) mcols(x[i1])$mpos == start(y[i2]) &
###       mcols(y[i2])$mpos == start(x[i1])
###
###   (D) isMateMinusStrand(x[i1]) == isMinusStrand(y[i2]) &
###       isMateMinusStrand(y[i2]) == isMinusStrand(x[i1])
###
###   (E) isFirstSegment(x[i1]) & isLastSegment(y[i2]) |
###       isFirstSegment(y[i2]) & isLastSegment(x[i1])
###
###   (F) isProperPair(x[i1]) == isProperPair(y[i2])
###
###   (G) isNotPrimaryRead(x[i1]) == isNotPrimaryRead(y[i2])

.checkMetadatacols <- function(arg, argname)
{
    if (!is(arg, "GappedAlignments"))
        stop("'", argname, "' must be a GappedAlignments object")
    if (is.null(names(arg)))
        stop("'", argname, "' must have names")
    arg_mcols <- mcols(arg)
    REQUIRED_COLNAMES <- c("flag", "mrnm", "mpos")
    if (!all(REQUIRED_COLNAMES %in% colnames(arg_mcols))) {
        colnames_in1string <-
            paste0("\"", REQUIRED_COLNAMES, "\"", collapse=", ")
        stop("required columns in 'mcols(", argname, ")': ",
             colnames_in1string)
    }
    arg_mcols
}

### 'names', 'flagbits', 'mrnm', and 'mpos', must all come from the same
### GappedAlignments object x.
### 'names': names(x).
### 'flagbits': integer matrix (of 0's and 1's) obtained with
###     bamFlagAsBitMatrix(mcols(x)$flag, bitnames=.MATING_FLAG_BITNAMES)
### 'mrnm': character vector or factor obtained with mcols(x)$mrnm
### 'mpos': integer vector obtained with mcols(x)$mpos
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

### Should return the same as:
###   args <- as.list(setNames(rep(TRUE, length(bitnames)), bitnames))
###   tmp <- do.call(scanBamFlag, args)
###   tmp[[2L]] - tmp[[1L]]
.makeFlagBitmask <- function(bitnames)
{
    bitpos <- match(bitnames, FLAG_BITNAMES)
    sum(as.integer(2L ^ (bitpos-1L)))
}

### All arguments must be atomic vectors of the same length N.
### The arguments prefixed with 'x_' describe a vector 'x' of N alignments.
### The arguments prefixed with 'y_' describe a vector 'y' of N alignments.
### Performs element-wise comparison of the N alignments in 'x' with the N
### alignments in 'y', and returns a logical vector of length N indicating
### whether conditions (B), (C), (D), (E), (F), and (G), are satisfied.
.isValidHit <- function(x_seqnames, x_start, x_mrnm, x_mpos, x_flag,
                        y_seqnames, y_start, y_mrnm, y_mpos, y_flag)
{
    D1_bitmask <- .makeFlagBitmask("isMateMinusStrand")
    D2_bitmask <- .makeFlagBitmask("isMinusStrand")
    E_bitmask <- .makeFlagBitmask("isFirstMateRead")
    FG_bitmask <- .makeFlagBitmask(c("isProperPair", "isNotPrimaryRead"))
    ## (B)
    x_mrnm == y_seqnames & y_mrnm == x_seqnames &
    ## (C)
      x_mpos == y_start & y_mpos == x_start &
    ## (D)
      (bitAnd(x_flag, D1_bitmask) != 0L) == (bitAnd(y_flag, D2_bitmask) != 0L) &
      (bitAnd(y_flag, D1_bitmask) != 0L) == (bitAnd(x_flag, D2_bitmask) != 0L) &
    ## (E)
      bitAnd(x_flag, E_bitmask) != bitAnd(y_flag, E_bitmask) &
    ## (F) & (G)
      bitAnd(x_flag, FG_bitmask) == bitAnd(y_flag, FG_bitmask)
}

### 3 equivalent implementations for this:
###   (a) x %in% x[duplicated(x)]
###   (b) duplicated(x) | duplicated(x, fromLast=TRUE)
###   (c) xx <- match(x, x); ans <- xx != seq_along(xx); ans[xx] <- ans; ans
### Comparing the 3 implementations on an integer vector of length 12 millions:
###   (a) is the most memory efficient;
###   (b) is a little bit faster than (a) (by only 8%) but uses between 12-14%
###       more memory;
###   (c) is as fast as (a) but uses about 30% more memory.
.hasDuplicates <- function(x)
{
    x %in% x[duplicated(x)]
}

.makeMateIdx <- function(x_hits, y_hits, x_len)
{
    oo <- IRanges:::orderInteger(y_hits, decreasing=TRUE)
    ans <- rep.int(NA_integer_, x_len)
    ans[x_hits[oo]] <- y_hits[oo]
    ans
}

### 'x_hits' and 'y_hits' must be 2 integer vectors of the same length N
### representing the N edges of a bipartite graph between the [1, x_len] and
### [1, y_len] intervals (the i-th edge being represented by (x[i], y[i])).
### Returns an integer vector F of length 'x_len' where F[k] is defined by:
###   - If there is no occurence of k in 'x', then F[k] = NA.
###   - If there is more than 1 occurence of k in 'x', then F[k] = 0.
###   - If there is exactly 1 occurence of k in 'x', at index i_k, then
###     F[k] = y[i_k].
### In addition, if more than valule of index k is associated to F[k], then
### F[k] is replaced by -F[k].
.makeMateIdx2 <- function(x_hits, y_hits, x_len)
{
    idx1 <- which(.hasDuplicates(y_hits))
    y_hits[idx1] <- - y_hits[idx1]
    idx2 <- which(.hasDuplicates(x_hits))
    y_hits[idx2] <- 0L
    ans <- rep.int(NA_integer_, x_len)
    ans[x_hits] <- y_hits
    ans
}

### Assumes that GappedAlignments object 'x':
###   (1) has clean group names i.e. .makeGappedAlignmentsGNames() would not
###       inject any NA in its names,
###   (2) is already ordered by group names (which are the same as its names).
### 'group.sizes' must be the same as 'runLength(Rle(names(x)))'.
.findMateAlignmentInChunk <- function(group.sizes,
                                      x_seqnames, x_start,
                                      x_mrnm, x_mpos, x_flag)
{
    hits <- IRanges:::makeAllGroupInnerHits(group.sizes, hit.type=1L)
    x_hits <- queryHits(hits)
    y_hits <- subjectHits(hits)
    valid_hits <- .isValidHit(x_seqnames[x_hits], x_start[x_hits],
                              x_mrnm[x_hits], x_mpos[x_hits], x_flag[x_hits],
                              x_seqnames[y_hits], x_start[y_hits],
                              x_mrnm[y_hits], x_mpos[y_hits], x_flag[y_hits])
    x_hits <- x_hits[valid_hits]
    y_hits <- y_hits[valid_hits]

    tmp <- x_hits
    x_hits <- c(x_hits, y_hits)
    y_hits <- c(y_hits, tmp)
    .makeMateIdx2(x_hits, y_hits, length(x_start))
}

.showGappedAlignmentsEltsWithMoreThan1Mate <- function(x, idx)
{
    if (length(idx) == 0L)
        return()
    cat("\n!! Found more than 1 mate for the following elements in 'x': ",
        paste(idx, collapse=", "),
        ".\n!! Details:\n!! ", sep="")
    GenomicRanges:::showGappedAlignments(x[idx],
                                         margin="!! ",
                                         with.classinfo=TRUE,
                                         print.seqlengths=FALSE)
    cat("!! ==> won't assign a mate to them!\n")
}

.dump_envir <- new.env(hash=TRUE, parent=emptyenv())
.dumpEnvir <- function() .dump_envir

flushDumpedAlignments <- function()
{
    objnames <- ls(envir=.dumpEnvir())
    rm(list=objnames, envir=.dumpEnvir())
}

.dumpAlignments <- function(x, idx)
{
    objnames <- ls(envir=.dumpEnvir())
    nobj <- length(objnames)
    if (nobj == 0L) {
        new_objname <- 1L
    } else {
        new_objname <- as.integer(objnames[nobj]) + 1L
    }
    new_objname <- sprintf("%08d", new_objname)
    assign(new_objname, x[idx], envir=.dumpEnvir())
}

countDumpedAlignments <- function()
{
    sum(unlist(eapply(.dumpEnvir(), length, USE.NAMES=FALSE)))
}

getDumpedAlignments <- function()
{
    objnames <- ls(envir=.dumpEnvir())
    args <- unname(mget(objnames, envir=.dumpEnvir()))
    do.call(c, args)
}

### Takes about 2.3 s and 170MB of RAM to mate 1 million alignments,
### and about 13 s and 909MB of RAM to mate 5 million alignments.
### So it's a little bit faster and more memory efficient than
### findMateAlignment2().
findMateAlignment <- function(x, verbose=FALSE)
{
    if (!isTRUEorFALSE(verbose))
        stop("'verbose' must be TRUE or FALSE")
    x_names <- names(x)
    if (is.null(x_names))
        stop("'x' must have names")
    x_mcols <- .checkMetadatacols(x, "x")
    ## flushDumpedAlignments() must be placed *after* the first reference to
    ## 'x', otherwise, when doing 'findMateAlignment(getDumpedAlignments())',
    ## the flushing would happen before 'x' is evaluated, causing 'x' to be
    ## evaluated to NULL.
    flushDumpedAlignments()
    x_flag <- x_mcols$flag
    bitnames <- c(.MATING_FLAG_BITNAMES, "isMinusStrand", "isMateMinusStrand")
    x_flagbits <- bamFlagAsBitMatrix(x_flag, bitnames=bitnames)
    x_mrnm <- as.character(x_mcols$mrnm)
    x_mpos <- x_mcols$mpos
    x_gnames <- .makeGappedAlignmentsGNames(x_names, x_flagbits,
                                            x_mrnm, x_mpos)
    x_seqnames <- as.character(seqnames(x))
    x_start <- start(x)

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
        chunk.x_seqnames <- x_seqnames[chunk.idx]
        chunk.x_start <- x_start[chunk.idx]
        chunk.x_mrnm <- x_mrnm[chunk.idx]
        chunk.x_mpos <- x_mpos[chunk.idx]
        chunk.x_flag <- x_flag[chunk.idx]
        if (verbose)
            message("Finding mates in chunk of ", chunk.length,
                    " alignments ... ", appendLF=FALSE)
        chunk.ans <- .findMateAlignmentInChunk(chunk.GS,
                                               chunk.x_seqnames,
                                               chunk.x_start,
                                               chunk.x_mrnm,
                                               chunk.x_mpos,
                                               chunk.x_flag)
        dumpme_idx <- which(chunk.ans <= 0L)
        if (length(dumpme_idx) != 0L) {
            .dumpAlignments(x, chunk.idx[dumpme_idx])
            chunk.ans[dumpme_idx] <- NA_integer_
        }
        if (verbose)
            message("OK")
        ans[chunk.idx] <- chunk.idx[chunk.ans]
        chunk.GIDX <- chunk.GIDX + NGROUP_BY_CHUNK
        chunk.offset <- chunk.offset + chunk.length
    }
    dump_count <- countDumpedAlignments()
    if (dump_count != 0L)
        warning("  ", dump_count, " alignments with ambiguous pairing ",
                "were dumped.\n    Use 'getDumpedAlignments()' to retrieve ",
                "them from the dump environment.")
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

### Takes about 2.8 s and 196MB of RAM to mate 1 million alignments,
### and about 19 s and 1754MB of RAM to mate 5 million alignments.
findMateAlignment2 <- function(x, y=NULL)
{
    x_names <- names(x)
    if (is.null(x_names))
        stop("'x' must have names")
    x_mcols <- .checkMetadatacols(x, "x")
    x_seqnames <- as.character(seqnames(x))
    x_start <- start(x)
    x_mrnm <- as.character(x_mcols$mrnm)
    x_mpos <- x_mcols$mpos
    x_flag <- x_mcols$flag
    bitnames <- c(.MATING_FLAG_BITNAMES, "isMinusStrand", "isMateMinusStrand")
    x_flagbits <- bamFlagAsBitMatrix(x_flag, bitnames=bitnames)
    x_gnames <- .makeGappedAlignmentsGNames(x_names, x_flagbits,
                                            x_mrnm, x_mpos)

    if (is.null(y)) {
        y_seqnames <- x_seqnames
        y_start <- x_start
        y_mrnm <- x_mrnm
        y_mpos <- x_mpos
        y_flag <- x_flag

        hits <- .findSelfMatches.character(x_gnames)
    } else {
        y_names <- names(y)
        if (is.null(y_names))
            stop("'y' must have names")
        y_mcols <- .checkMetadatacols(y, "y")
        y_seqnames <- as.character(seqnames(y))
        y_start <- start(y)
        y_mrnm <- as.character(y_mcols$mrnm)
        y_mpos <- y_mcols$mpos
        y_flag <- y_mcols$flag
        y_flagbits <- bamFlagAsBitMatrix(y_flag, bitnames=bitnames)
        y_gnames <- .makeGappedAlignmentsGNames(y_names, y_flagbits,
                                                y_mrnm, y_mpos)

        hits <- .findMatches(x_gnames, y_gnames, incomparables=NA_character_)
    }

    x_hits <- queryHits(hits)
    y_hits <- subjectHits(hits)
    valid_hits <- .isValidHit(x_seqnames[x_hits], x_start[x_hits],
                              x_mrnm[x_hits], x_mpos[x_hits], x_flag[x_hits],
                              y_seqnames[y_hits], y_start[y_hits],
                              y_mrnm[y_hits], y_mpos[y_hits], y_flag[y_hits])
    x_hits <- x_hits[valid_hits]
    y_hits <- y_hits[valid_hits]

    if (is.null(y)) {
        tmp <- x_hits
        x_hits <- c(x_hits, y_hits)
        y_hits <- c(y_hits, tmp)
    }
    ans <- .makeMateIdx2(x_hits, y_hits, length(x))
    if (any(ans <= 0L, na.rm=TRUE)) {
        more_than_1_mate_idx <- which(ans == 0L)
        .showGappedAlignmentsEltsWithMoreThan1Mate(x, more_than_1_mate_idx)
        ans[ans <= 0L] <- NA_integer_
    }
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
    .isFirstSegment.integer(mcols(x)$flag)

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
    .isLastSegment.integer(mcols(x)$flag)

### 'x' must be a GappedAlignments objects.
makeGappedAlignmentPairs <- function(x, use.names=FALSE, use.mcols=FALSE)
{
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    if (!isTRUEorFALSE(use.mcols)) {
        if (!is.character(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE or a character vector ",
                 "specifying the metadata columns to propagate")
        if (!all(use.mcols %in% colnames(mcols(x))))
            stop("'use.mcols' must be a subset of 'colnames(mcols(x))'")
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
    x_is_proper <- as.logical(bamFlagAsBitMatrix(mcols(x)$flag,
                                                 bitnames="isProperPair"))
    ans_is_proper <- x_is_proper[first_idx]

    ## Drop pairs with discordant seqnames or strand.
    idx_is_discordant <- (as.character(seqnames(x)[first_idx]) !=
                          as.character(seqnames(x)[last_idx])) |
                         (as.character(strand(x)[first_idx]) ==
                          as.character(strand(x)[last_idx]))
    if (any(idx_is_discordant) != 0L) {
        nb_discordant_proper <- sum(ans_is_proper[idx_is_discordant])
        if (nb_discordant_proper != 0L) {
            ratio <- 100.0 * nb_discordant_proper / sum(idx_is_discordant)
            warning(ratio, "% of the pairs with discordant seqnames or ",
                    "strand were flagged\n",
                    "  as proper pairs by the aligner. Dropping them anyway.")
        }
        keep <- -which(idx_is_discordant)
        first_idx <- first_idx[keep]
        last_idx <- last_idx[keep]
        ans_is_proper <- ans_is_proper[keep]
    }

    ## The big split!
    ans_first <- x[first_idx]
    ans_last <- x[last_idx]
    ans_names <- NULL
    if (use.names)
        ans_names <- names(ans_first)
    names(ans_first) <- names(ans_last) <- NULL
    if (is.character(use.mcols)) {
        mcols(ans_first) <- mcols(ans_first)[use.mcols]
        mcols(ans_last) <- mcols(ans_last)[use.mcols]
    } else if (!use.mcols) {
        mcols(ans_first) <- mcols(ans_last) <- NULL
    }
    GappedAlignmentPairs(ans_first, ans_last, ans_is_proper, names=ans_names)
}

