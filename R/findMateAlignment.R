### =========================================================================
### findMateAlignment()
### -------------------------------------------------------------------------
### Find the mate of each element in a GappedAlignments object 'x'.
### 'x[i1]' and 'x[i2]' are considered mate iff:
###   (A) names(x[i1]) == names(x[i2])
###   (B) elementMetadata(x[i1])$mrnm == seqnames(x[i2])
###   (C) elementMetadata(x[i1])$mpos == start(x[i2])
###   (D) isFirstSegment(x[i1]) & isLastSegment(x[i2]) |
###       isLastSegment(x[i1]) & isFirstSegment(x[i2])
### We do NOT look at the "isize" column (TLEN field in SAM Spec).

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
        subject_low2high <- IRanges:::.makeLow2highFromHigh2low(
                                high2low(subject))
        subject_hits0 <- m0[query_hits0]
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

### TODO: Maybe put this in IRanges. Maybe operate on a Hits object?
.hasReverseHit <- function(query_hits, subject_hits)
{
    if (!is.integer(query_hits) || !is.integer(subject_hits))
        stop("'query_hits' and 'subject_hits' must be integer vectors")
    if (length(query_hits) != length(subject_hits))
        stop("'query_hits' and 'subject_hits' must have the same length")
    is_dup <- IRanges:::duplicatedIntegerPairs(c(subject_hits, query_hits),
                                               c(query_hits, subject_hits))
    is_dup[length(query_hits)+seq_len(length(query_hits))]
}

### Takes about 10 sec and 250MB of RAM to mate 1 million alignments.
findMateAlignment <- function(x)
{
    if (!is(x, "GappedAlignments"))
        stop("'x' must be a GappedAlignments object")
    x_eltmetadata <- elementMetadata(x)
    X_REQUIRED_COLS <- c("flag", "mrnm", "mpos")
    if (!all(X_REQUIRED_COLS %in% colnames(x_eltmetadata))) {
        cols_in_1string <- paste("\"", X_REQUIRED_COLS, "\"", sep="",
                                 collapse=", ")
        stop("required columns in 'elementMetadata(x)': ", cols_in_1string)
    }

    x_names <- names(x)
    ## Before we pass 'x_names' to .findMatches() we inject NAs at positions
    ## corresponding to alignments with an NA 'x_mpos' (i.e. PNEXT = 0), or to
    ## unpaired reads, or to reads that are neither first or last mate.
    x_mpos <- x_eltmetadata$mpos
    x_flag <- x_eltmetadata$flag
    x_flagbits <- bamFlagAsBitMatrix(x_flag)
    x_is_first <- .isFirstSegment.matrix(x_flagbits)
    x_is_last <- .isLastSegment.matrix(x_flagbits)
    x_ignore <- which(is.na(x_mpos) | !(x_is_first | x_is_last))
    x_names[x_ignore] <- NA_integer_

    hits <- .findMatches(x_names, x_names, incomparables=NA_character_)
    x_hits <- queryHits(hits)
    y_hits <- subjectHits(hits)

    ## Keep hits that satisfy condition (C).
    y_start <- start(x)
    hit_is_C <- x_mpos[x_hits] == y_start[y_hits]
    x_hits <- x_hits[hit_is_C]
    y_hits <- y_hits[hit_is_C]

    ## Keep hits that sattisfy condition (B).
    x_mrnm <- as.character(x_eltmetadata$mrnm)
    y_seqnames <- as.character(seqnames(x))
    hit_is_B <- x_mrnm[x_hits] == y_seqnames[y_hits]
    x_hits <- x_hits[hit_is_B]
    y_hits <- y_hits[hit_is_B]

    ## Keep hits that satisfy condition (D).
    hit_is_D <- x_is_first[x_hits] & x_is_last[y_hits] |
                x_is_last[x_hits] & x_is_first[y_hits]
    x_hits <- x_hits[hit_is_D]
    y_hits <- y_hits[hit_is_D]

    ## Keep hits that have a reverse hit.
    has_rev_hit <- .hasReverseHit(x_hits, y_hits)
    x_hits <- x_hits[has_rev_hit]
    y_hits <- y_hits[has_rev_hit]

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

