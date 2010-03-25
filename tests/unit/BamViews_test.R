.BamViews_ok <-
    function(v, dim=c(0,0), bamRanges=GRanges(),
             bamSamples=DataFrame(row.names=character()),
             bamExperiment=list())
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
    rd <- GRanges(rep(c("chr1", "chr2"), each=5),
                  IRanges(c(1:5, 101:105), c(11:15, 111:115)),
                  Values=rev(seq_len(ni)))

    sd0 <- new("DataFrame", nrows=nj,
               rownames=make.unique(basename(fls)))
    sd1 <-DataFrame(Score=seq_len(nj),
                     row.names=make.unique(basename(fls)))

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
    rl2 <- GRanges(rep(c("chr1", "chr2"), each=5),
                   IRanges(c(1:5, 101:105), c(11:15, 111:115)))
    nj1 <- 4L
    fls1 <- rep("foo", nj1)
    sd1 <- DataFrame(Value=rev(seq_len(nj1)))
    ni2 <- length(rl2)

    rd2 <- rl2
    elementMetadata(rd2)[["Count"]] <- rev(seq_len(ni2))
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

test_BamViews_bamIndicies <- function()
{
    bv <- BamViews()
    checkIdentical(character(0), bamIndicies(bv))

    ## copy fls as bamIndicies
    fls <- c(tempfile(), tempfile())    # does not exist
    bv <- BamViews(fls)
    checkIdentical(fls, bamPaths(bv))
    checkIdentical(fls, bamIndicies(bv))

    ## keep bamPaths, bamIndicies separate
    ifls <- c(tempfile(), tempfile())    # does not exist
    bv <- BamViews(fls, ifls)
    checkIdentical(fls, bamPaths(bv))
    checkIdentical(ifls, bamIndicies(bv))

    ## subsetting
    checkIdentical(fls[2], bamPaths(bv[,2]))
    checkIdentical(ifls[1], bamIndicies(bv[,1]))
}

test_BamViews_auto.range <- function()
{
    fl <- tempfile()                    # does not exist
    checkIdentical(GRanges(), bamRanges(BamViews(fl)))

    bv <- msg <- NULL
    suppressWarnings({
        withCallingHandlers({
            bv <- BamViews(fl, auto.range=TRUE)
        }, warning=function(w) {
            msg <<- conditionMessage(w)
        })
    })
    checkIdentical(GRanges(), bamRanges(bv))
    checkIdentical("some files do not exist; bamRanges not defined",
                   msg)

    fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
    bv <- BamViews(c(fl, fl), auto.range=FALSE)
    checkIdentical(GRanges(), bamRanges(bv))

    bv <- BamViews(c(fl, fl), auto.range=TRUE)
    checkIdentical(GRanges(c("seq1", "seq2"), IRanges(1, c(1575, 1584))),
                   bamRanges(bv))
}
    
test_BamViews_readBamGappedAlignments <- function()
{
    checkTrue(validObject(readBamGappedAlignments(BamViews())))

    fl <- c(system.file("extdata", "ex1.bam", package="Rsamtools"),
            file.path("cases", "ex1_shuf1000.bam"))
    bv <- BamViews(fl, auto.range=TRUE)
    rng <- bamRanges(bv)
    aln <- readBamGappedAlignments(bv)
    checkEquals(length(bamPaths(bv)), length(aln))

    fl <- c(fl, tempfile())
    bv <- BamViews(fl, bamRanges=rng)
    msg <- NULL
    suppressWarnings({
        tryCatch({
            aln <- readBamGappedAlignments(bv)
        }, error=function(err) {
            msg <<- conditionMessage(err)
        })
    })
    tst <- sprintf("'readBamGappedAlignments' failed on '%s'",
                   names(bv)[3])
    checkIdentical(tst, msg)
}
