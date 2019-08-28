fl <- system.file("extdata", "ex1.cram", package = "Rsamtools")
idx <- paste0(fl, ".crai")

test_BamFile <- function()
{
    bfl <- BamFile(fl)
    checkIdentical(fl, path(bfl))
    checkIdentical(idx, index(bfl))
    checkTrue(!isOpen(bfl))

    bfl <- BamFile(fl, idx)
    checkIdentical(fl, path(bfl))
    checkIdentical(idx, index(bfl))
    checkTrue(!isOpen(bfl))

    open(bfl)
    checkTrue(isOpen(bfl))

    close(bfl)
    checkTrue(!isOpen(bfl))


    bfl <- BamFile(tmpfl <- tempfile(fileext=".crai"))
    checkIdentical(tmpfl, path(bfl))
    checkIdentical(NA_character_, index(bfl))
    checkException(open(bfl), silent = TRUE)
}

test_seqinfo <- function()
{
    bfl <- BamFile(fl)
    checkTrue(validObject(si <- seqinfo(bfl)))
    checkIdentical("21", seqnames(si))
    checkIdentical(c(`21` = 46944323L), seqlengths(si))
    checkIdentical(c(`21` = NA), isCircular(si))
    checkIdentical(c(`21` = NA_character_), genome(si))
}

test_idxstatsBam <- function()
{
    target <- structure(
        list(
            seqnames = structure(1:2, .Label = c("21", "*"), class = "factor"),
            seqlength = c(46944323L, 0L),
            mapped = c(1673, 0),
            unmapped = c(22, 0)),
        class = "data.frame",
        row.names = c(NA, -2L)
    )
Bma    checkIdentical(target, idxstatsBam(fl, idx))
    checkIdentical(target, idxstatsBam(fl))

    ## idxstats changes data read from cram; does it reset correctly?
    bf <- BamFile(fl)
    open(bf)
    checkIdentical(target, idxstatsBam(bf))
    checkIdentical(69957, countBam(bf)$nucleotides)
}

test_countBam <- function()
{
    ## entire file
    target <- data.frame(
        space = NA,
        start = NA, end = NA, width = NA,
        file = "ex1.cram",
        records = 1695L, nucleotides = 69957
    )
    checkIdentical(target, countBam(fl))

    bfl <- BamFile(fl)
    checkIdentical(target, countBam(bfl))

    open(bfl)
    checkIdentical(target, countBam(bfl))
    close(bfl)

    bfl <- BamFile(fl, idx)
    checkIdentical(target, countBam(bfl))

    ## entire sequence
    param <- ScanBamParam(which = GRanges("21:1-46944323"))
    target <- data.frame(
        space = "21",
        start = 1L, end = 46944323L, width = 46944323L,
        file = "ex1.cram",
        records = 1695L, nucleotides = 69957
    )
    checkIdentical(target, countBam(bfl, param = param))

    ## specific coordinates
    start <- 33049956L
    end <- 33052045L
    param <- ScanBamParam(which = GRanges("21", IRanges(start, end)))
    target <- data.frame(
        space = "21",
        start = start, end = end, width = end - start + 1L,
        file = "ex1.cram",
        records = 1695L, nucleotides = 69957
    )
    checkIdentical(target, countBam(bfl, param = param))
    
    ## duplicate ranges
    param <- ScanBamParam(which = GRanges("21", IRanges(start, end))[c(1, 1)])
    target1 <- target[c(1, 1),]
    rownames(target1) <- NULL
    checkIdentical(target1, countBam(bfl, param = param))

    ## overlapping ranges
    mid <- as.integer(start + (end - start) / 2)
    param <- ScanBamParam(
        which = GRanges("21", IRanges(c(start, mid), end))
    )
    target <- data.frame(
        space = "21",
        start = c(start, mid), end = end,
        width = c(end - start + 1L, end - mid + 1L),
        file = "ex1.cram",
        records = c(1695L, 863L),
        nucleotides = c(69957, 35300)
    )
    checkIdentical(target, countBam(bfl, param = param))

    param <- ScanBamParam(
        which = GRanges("21", IRanges(c(start, mid), end))[2:1]
    )
    target1 <- target[2:1,]
    rownames(target1) <- NULL
    checkIdentical(target1, countBam(bfl, param = param))
}

test_scanBam <- function() {
    res <- scanBam(fl)
    checkIdentical(13L, lengths(res))
    target <- c(
        qname = 1695L, flag = 1695L, rname = 1695L, strand = 1695L, 
        pos = 1695L, qwidth = 1695L, mapq = 1695L, cigar = 1695L, mrnm = 1695L, 
        mpos = 1695L, isize = 1695L, seq = 1695L, qual = 1695L
    )
    checkIdentical(target, lengths(res[[1]]))
    checkIdentical(res, scanBam(BamFile(fl)))
    checkIdentical(res, scanBam(BamFile(fl, idx)))

    start <- 33049956L
    end <- 33052045L
    mid <- floor(start + (end - start) / 2)
    param <- ScanBamParam(
        which = GRanges("21", IRanges(start, mid)),
        what = scanBamWhat()
    )
    res <- scanBam(BamFile(fl), param = param)
    checkIdentical(c(`21:33049956-33051000` = 13L), lengths(res))
    target <- c(
        qname = 864L, flag = 864L, rname = 864L, strand = 864L, pos = 864L, 
        qwidth = 864L, mapq = 864L, cigar = 864L, mrnm = 864L, mpos = 864L, 
        isize = 864L, seq = 864L, qual = 864L
    )
    checkIdentical(target, lengths(res[[1]]))
}
