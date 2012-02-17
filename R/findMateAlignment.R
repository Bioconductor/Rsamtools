### =========================================================================
### findMateAlignment()
### -------------------------------------------------------------------------
### For each element in GappedAlignments object 'x', finds its mate in
### GappedAlignments object 'y'.
### 'x[i1]' and 'y[i2]' are considered mate iff:
###   (A) names(x[i1]) == names(y[i2])
###   (B) elementMetadata(x[i1])$mrnm == seqnames(y[i2]) &
###       elementMetadata(y[i2])$mrnm == seqnames(x[i1])
###   (C) elementMetadata(x[i1])$mpos == start(y[i2]) &
###       elementMetadata(y[i2])$mpos == start(x[i1])
###   (D) isFirstSegment(x[i1]) & isLastSegment(y[i2]) |
###       isFirstSegment(y[i2]) & isLastSegment(x[i1])
### We do NOT look at the "isize" column (TLEN field in SAM Spec).

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
    .isFirstSegment.matrix(bamFlagAsBitMatrix(flag))

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
    .isLastSegment.matrix(bamFlagAsBitMatrix(flag))

.isLastSegment.GappedAlignments <- function(x)
    .isLastSegment.integer(elementMetadata(x)$flag)

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

### Use to find all self matches in 'x'. A little bit faster than
### 'findMatches(x, x, incomparables=NA_character_)' but not so much.
.findMatches2.character <- function(x)
{
    xx <- match(x, x, nomatch=0L, incomparables=NA_character_)
    xxo <- IRanges:::orderInteger(xx)
    xx2 <- Rle(xx[xxo])
    shift0 <- 0L
    GS <- runLength(xx2)
    if (runValue(xx2)[1L] == 0L) {
        shift0 <- GS[1L]
        GS <- GS[-1L]
    }
    ans <- IRanges:::makeAllGroupInnerHits(GS)
    query_hits <- xxo[ans@queryHits + shift0]
    subject_hits <- xxo[ans@subjectHits + shift0]
    oo <- IRanges:::orderIntegerPairs(query_hits, subject_hits)
    ans@queryHits <- query_hits[oo]
    ans@subjectHits <- subject_hits[oo]
    ans@queryLength <- ans@subjectLength <- length(x)
    ans
}

### Takes about 11 sec and 260MB of RAM to mate 1 million alignments.
findMateAlignment <- function(x, y=NULL)
{
    x_eltmetadata <- .checkElementMetadata(x, "x")
    x_names <- names(x)
    x_seqnames <- as.character(seqnames(x))
    x_start <- start(x)
    x_flag <- x_eltmetadata$flag
    x_mrnm <- as.character(x_eltmetadata$mrnm)
    x_mpos <- x_eltmetadata$mpos

    ## Before we pass 'x_names' to .findMatches() we inject NAs at positions
    ## corresponding to alignments with an NA in 'x_mrnm' or 'x_mpos' (i.e.
    ## RNEXT = '*' or PNEXT = 0), or to unpaired reads, or to reads that are
    ## neither first or last mate.
    x_flagbits <- bamFlagAsBitMatrix(x_flag)
    x_is_first <- .isFirstSegment.matrix(x_flagbits)
    x_is_last <- .isLastSegment.matrix(x_flagbits)
    x_ignore <- which(is.na(x_mrnm) |
                      is.na(x_mpos) |
                      !(x_is_first | x_is_last))
    x_names[x_ignore] <- NA_integer_

    if (is.null(y)) {
        y_seqnames <- x_seqnames
        y_start <- x_start
        y_mrnm <- x_mrnm
        y_mpos <- x_mpos
        y_is_first <- x_is_first
        y_is_last <- x_is_last
        hits <- .findMatches(x_names, x_names, incomparables=NA_character_)
        #hits <- .findMatches2.character(x_names)
    } else {
        y_eltmetadata <- .checkElementMetadata(y, "y")
        y_names <- names(y)
        y_seqnames <- as.character(seqnames(y))
        y_start <- start(y)
        y_flag <- y_eltmetadata$flag
        y_mrnm <- as.character(y_eltmetadata$mrnm)
        y_mpos <- y_eltmetadata$mpos

        ## Before we pass 'y_names' to .findMatches() we inject NAs at positions
        ## corresponding to alignments with an NA in 'y_mrnm' or 'y_mpos' (i.e.
        ## RNEXT = '*' or PNEXT = 0), or to unpaired reads, or to reads that are
        ## neither first or last mate.
        y_flagbits <- bamFlagAsBitMatrix(y_flag)
        y_is_first <- .isFirstSegment.matrix(y_flagbits)
        y_is_last <- .isLastSegment.matrix(y_flagbits)
        y_ignore <- which(is.na(y_mrnm) |
                          is.na(y_mpos) |
                          !(y_is_first | y_is_last))
        y_names[y_ignore] <- NA_integer_
        hits <- .findMatches(x_names, y_names, incomparables=NA_character_)
    }

    x_hits <- queryHits(hits)
    y_hits <- subjectHits(hits)

    ## Keep hits that satisfy condition (C).
    hit_is_C <- x_mpos[x_hits] == y_start[y_hits] &
                y_mpos[y_hits] == x_start[x_hits]
    x_hits <- x_hits[hit_is_C]
    y_hits <- y_hits[hit_is_C]

    ## Keep hits that sattisfy condition (B).
    hit_is_B <- x_mrnm[x_hits] == y_seqnames[y_hits] &
                y_mrnm[y_hits] == x_seqnames[x_hits]
    x_hits <- x_hits[hit_is_B]
    y_hits <- y_hits[hit_is_B]

    ## Keep hits that satisfy condition (D).
    hit_is_D <- x_is_first[x_hits] & y_is_last[y_hits] |
                y_is_first[y_hits] & x_is_last[x_hits]
    x_hits <- x_hits[hit_is_D]
    y_hits <- y_hits[hit_is_D]

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

makeGappedAlignmentPairs <- function(x, use.names=FALSE)
{
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    mate <- findMateAlignment(x)
    is_first <- .isFirstSegment.GappedAlignments(x)
    is_last <- .isLastSegment.GappedAlignments(x)
    first_idx <- which(!is.na(mate) & is_first)
    last_idx <- mate[first_idx]

    ## Fundamental property of the 'mate' vector: it's a permutation of order
    ## 2 and with no fixed point on the set of indices for which 'mate' is
    ## not NA.
    ## Checking there are no fixed points:
    if (!all(first_idx != last_idx))
        stop("findMateAlignment() returned an invalid 'mate' vector")
    ## Checking order 2 (i.e. permuting a 2nd time brings back the original
    ## set of indices):
    if (!identical(mate[last_idx], first_idx))
        stop("findMateAlignment() returned an invalid 'mate' vector")
    ## One more sanity check:
    if (!all(is_last[last_idx]))
        stop("findMateAlignment() returned an invalid 'mate' vector")

    first <- x[first_idx]
    last <- x[last_idx]
    first_flagbits <- bamFlagAsBitMatrix(elementMetadata(first)$flag)
    first_is_proper <- as.logical(first_flagbits[ , "isProperPair"])
    last_flagbits <- bamFlagAsBitMatrix(elementMetadata(last)$flag)
    last_is_proper <- as.logical(last_flagbits[ , "isProperPair"])
    if (!identical(first_is_proper, last_is_proper))
        stop("for some pairs, the 2 mates don't have ",
             "the same \"isProperPair\" bit.\n",
             "  Maybe the BAM file they are coming from was corrupted or ",
             "generated by\n  an aligner that doesn't follow the SAM Spec?")
    ans_names <- NULL
    if (use.names)
        ans_names <- names(first)
    names(first) <- names(last) <- NULL
    GappedAlignmentPairs(first, last, first_is_proper, names=ans_names)
}

