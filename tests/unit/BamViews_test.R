.BamViews_ok <-
    function(v, dim=c(0,0), bamRanges=RangedData(),
             bamSamples=DataFrame(), bamExperiment=list())
{
    checkTrue(validObject(v))
    checkEquals(dim, dim(v))
    checkIdentical(bamRanges, bamRanges(v))
    checkIdentical(bamSamples, bamSamples(v))
    checkIdentical(bamExperiment, bamExperiment(v))
}

test_BamViews_constructors <- function()
{
    .BamViews_ok(BamViews(), c(0, 0))
    ni <- 10L; nj <- 5L
    fls <- rep("foo", nj)
    rd <- RangedData(RangesList(chr1=IRanges(1:5, 11:15),
                                chr2=IRanges(101:105, 111:115)),
                     Values=rev(seq_len(ni)))
    sd0 <- new("DataFrame", nrows=nj)
    sd1 <- DataFrame(Score=seq_len(nj))

    .BamViews_ok(BamViews(fls), dim=c(0, nj), bamSamples=sd0)
    .BamViews_ok(BamViews(fls, bamRanges=rd), dim=c(ni, nj),
                 bamRanges=rd, bamSamples=sd0)
    .BamViews_ok(BamViews(fls, bamRanges=rd, bamSamples=sd0),
                 dim=c(ni, nj), bamRanges=rd, bamSamples=sd0)
    .BamViews_ok(BamViews(fls, bamRanges=rd, bamSamples=sd1),
                 dim=c(ni, nj), bamRanges=rd, bamSamples=sd1)
}

test_BamViews_subset <- function()
{
    rl2 <- RangesList(chr1=IRanges(1:5, 11:15),
                      chr2=IRanges(101:105, 111:115))
    nj1 <- 4L
    fls1 <- rep("foo", nj1)
    sd1 <- DataFrame(Value=rev(seq_len(nj1)))
    ni2 <- sum(sapply(rl2, length))

    rd2 <- RangedData(rl2, Count=rev(seq_len(ni2)))
    ## rows
    v <- BamViews(bamPaths=fls1, bamRanges=rd2, bamSamples=sd1)
    .BamViews_ok(v, c(ni2, nj1), bamRanges=rd2, bamSamples=sd1)
    .BamViews_ok(v[TRUE,], c(ni2, nj1), bamRanges=rd2, bamSamples=sd1)
    i_idx <- c(FALSE, TRUE)
    .BamViews_ok(v[i_idx,], c(ni2/2, nj1), bamRanges=rd2[i_idx,],
                 bamSamples=sd1)
    i_idx <- c(2, 4)
    .BamViews_ok(v[i_idx,],c(length(i_idx), nj1), bamRanges=rd2[i_idx,],
                 bamSamples=sd1)
    ## columns
    .BamViews_ok(v[,TRUE], c(ni2, nj1), bamRanges=rd2, bamSamples=sd1)
    j_idx <- c(FALSE, TRUE, FALSE, TRUE)
    .BamViews_ok(v[,j_idx], c(ni2, nj1/2), bamRanges=rd2,
                 bamSamples=sd1[j_idx,,drop=FALSE])
    j_idx <- c(2, 4)
    .BamViews_ok(v[,j_idx], c(ni2, nj1/2), bamRanges=rd2,
                 bamSamples=sd1[j_idx,,drop=FALSE])
}

test_BamViews_subset_RangesList <- function()
{
    v <- BamViews()
    msg <- NULL
    checkException(tryCatch(v[RangesList(),], error=function(err) {
        msg <<- conditionMessage(err)
        stop(err)
    }), silent=TRUE)
    checkIdentical("'[,BamViews,RangesList-method' not yet supported", msg)
##     rl0 <- RangesList()
##     rl1 <- RangesList(chr1=IRanges(1:5, 11:15))
##     rl2 <- RangesList(chr1=IRanges(1:5, 11:15),
##                       chr2=IRanges(101:105, 111:115))
##     nj <- 4L
##     fls <- rep("foo", nj)
##     sd0 <- DataFrame()
##     sd1 <- DataFrame(Value=rev(seq_len(nj)))
##     ni0 <- 0L
##     ni1 <- sum(sapply(rl1, length))
##     ni2 <- sum(sapply(rl2, length))

##     v <- BamViews(bamPaths=fls, bamSamples=sd1)
##     .BamViews_ok(v[rl0,], dim=c(ni0, nj), bamRanges=RangedData(),
##                  bamSamples=sd1)
##     .BamViews_ok(v[rl1,], dim=c(ni1, nj), bamRanges=RangedData(rl1),
##                  bamSamples=sd1)
##     .BamViews_ok(v[rl2,], dim=c(ni2, nj), bamRanges=RangedData(rl2),
##                  bamSamples=sd1)

##     rd1 <- RangedData(rl1, Values=seq_len(ni1))
##     v <- BamViews(bamPaths=fls, bamSamples=sd1, bamRanges=rd1)
##     .BamViews_ok(v[rl0,], dim=c(ni1, nj), bamRanges=rd1,
##                  bamSamples=sd1)
##     .BamViews_ok(v[rl1,], dim=c(ni1 * ni1, nj),
##                  bamRanges=bamRanges(v)[rep(seq_len(ni1), each=ni1),],
##                  bamSamples=sd1p)
##     .BamViews_ok(v[rl2,], dim=c(ni1, nj), bamRanges=rd1,
##                  bamSamples=sd1)
}
