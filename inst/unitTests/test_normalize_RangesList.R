.normalizeRangesList <- Rsamtools:::.normalizeRangesList

test_normalizeRangesList <- function()
{
    rl <- RangesList()
    checkIdentical(rl, .normalizeRangesList(rl))

    rl <- RangesList(IRanges())
    checkEquals(RangesList(as(IRanges(), "NormalIRanges")),
                .normalizeRangesList(rl))

    rl <- RangesList(IRanges(1,2))
    checkEquals(RangesList(as(IRanges(1, 2), "NormalIRanges")),
                .normalizeRangesList(rl))

    rl <- RangesList(IRanges(1:2,11:12))
    checkEquals(RangesList(as(IRanges(1, 12), "NormalIRanges")),
                .normalizeRangesList(rl))

    rl <- RangesList(IRanges(), IRanges())
    checkEquals(RangesList(as(IRanges(), "NormalIRanges")),
                .normalizeRangesList(rl))

    rl <- RangesList(x=IRanges(), y=IRanges())
    rln <- endoapply(rl, as, "NormalIRanges")
    checkEquals(rln, .normalizeRangesList(rl))

    rl <- RangesList(x=IRanges(), x=IRanges())
    rln <- endoapply(rl[1], as, "NormalIRanges")
    checkEquals(rln, .normalizeRangesList(rl))

    rl <- RangesList(x=IRanges(), IRanges())
    rln <- endoapply(rl, as, "NormalIRanges")
    checkEquals(rln, .normalizeRangesList(rl))

    rl <- RangesList(x=IRanges(1:2, 11:12), y=IRanges(1:2, 11:12))
    rln <- endoapply(rl, as, "NormalIRanges")
    checkEquals(rln, .normalizeRangesList(rl))

    rl <- RangesList(x=IRanges(1:2, 11:12), x=IRanges(1:2, 11:12))
    rln <- RangesList(x=as(IRanges(1, 12), "NormalIRanges"))
    checkEquals(rln, .normalizeRangesList(rl))

    rl <- RangesList(x=IRanges(1:2, 11:12), x=IRanges(3:4, 13:14))
    rln <- RangesList(x=as(IRanges(1, 14), "NormalIRanges"))
    checkEquals(rln, .normalizeRangesList(rl))

    rl <- RangesList(x=IRanges(1:2, 11:12), x=IRanges(3:4, 13:14),
                     y=IRanges(1:2, 11:12))
    rln <- RangesList(x=as(IRanges(1, 14), "NormalIRanges"),
                      y=as(IRanges(1,12), "NormalIRanges"))
    checkEquals(rln, .normalizeRangesList(rl))
}
