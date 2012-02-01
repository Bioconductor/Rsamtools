fl <- system.file("extdata", "ce2dict1.fa", package="Rsamtools")

test_FaFile_openclose <- function()
{
    fa <- FaFile(fl)
    checkIdentical(FALSE, isOpen(fa))
    checkIdentical(TRUE, isOpen(open(fa)))
    checkIdentical(FALSE, isOpen(close(fa)))
}

test_FaFile_emptyfile <- function()
{
    fl <- tempfile()
    file.create(fl)
    checkException(suppressWarnings(open(fa <- FaFile(fl))),
                   silent=TRUE)
}
    
test_FaFile_emptyid <- function()
{
    fl <- tempfile()
    cat(">\nACTA", file=fl)
    open(fa <- FaFile(fl))
    close(fa)
}

test_FaFile_scanFaIndex <- function()
{
    .checkIdx <- function(idx) {
        checkTrue(is(idx, "GRanges"))
        checkIdentical(5L, length(idx))
        checkIdentical(116L, sum(width(idx)))
    }
    fa <- FaFile(fl)
    .checkIdx(scanFaIndex(fa))
    .checkIdx(scanFaIndex(fl))
}

test_FaFile_count <- function()
{
    checkIdentical(5L, countFa(open(FaFile(fl))))
    checkIdentical(5L, countFa(fl))
}

test_FaFile_scanFa <- function()
{
    .checkRes <- function(res) {
        checkTrue(is(res, "DNAStringSet"))
        checkIdentical(5L, length(idx))
        checkIdentical(116L, sum(width(idx)))
    }
    fa <- open(FaFile(fl))
    idx <- scanFaIndex(fa)[1:5]
    .checkRes(scanFa(fa, idx))
    .checkRes(scanFa(fl, idx))

    ## scanFa,*,missing-methods
    checkTrue(validObject(scanFa(fa)))
    checkTrue(validObject(scanFa(fl)))

    ## GRanges
    exp0 <- subseq(scanFa(fa)[c(1, 3)], 5, 9)
    gr <- GRanges(c("pattern01", "pattern03"), IRanges(5, width=5))
    checkIdentical(as.character(exp0), as.character(scanFa(fa, gr)))
    checkIdentical(as.character(DNAStringSet()),
                   unname(as.character(scanFa(fa, GRanges()))))

    ## scanFa ignores strand
    strand(gr) <- c("-", "+")
    checkIdentical(as.character(exp0), as.character(scanFa(fa, gr)))

    ## RangedData, RangesList
    rd <- as(gr, "RangedData")
    checkIdentical(as.character(exp0), as.character(scanFa(fa, rd)))
    rl <- as(gr, "RangesList")
    checkIdentical(as.character(exp0), as.character(scanFa(fa, rd)))
}

test_FaFile_getSeq <- function()
{
    fa <- open(FaFile(fl))
    exp0 <- subseq(scanFa(fa)[c(1, 3)], 5, 9)

    gr <- GRanges(c("pattern01", "pattern03"), IRanges(5, width=5))
    checkIdentical(as.character(exp0), as.character(getSeq(fa, gr)))

    ## '-' strand is reverse complement
    strand(gr) <- c("-", "+")
    exp <- exp0
    exp[1] <- reverseComplement(exp[1])
    checkIdentical(as.character(exp), as.character(getSeq(fa, gr)))

    ## RangedData, RangesList
    rd <- as(gr, "RangedData")
    checkIdentical(as.character(exp0), as.character(getSeq(fa, rd)))
    rl <- as(gr, "RangesList")
    checkIdentical(as.character(exp0), as.character(getSeq(fa, rd)))
}
