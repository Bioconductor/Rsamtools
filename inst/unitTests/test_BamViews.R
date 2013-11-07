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

    sd0 <- DataFrame(row.names=make.unique(basename(fls)))
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
    mcols(rd2)[["Count"]] <- rev(seq_len(ni2))
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

test_BamViews_bamIndicies <- function()
{
    bv <- BamViews()
    checkIdentical(setNames(character(0), character(0)), bamIndicies(bv))

    ## copy fls as bamIndicies
    fls <- c(tempfile(), tempfile())    # does not exist
    bv <- BamViews(fls)
    fls <- setNames(fls, names(bv))
    checkIdentical(fls, bamPaths(bv))
    checkIdentical(fls, bamIndicies(bv))

    ## keep bamPaths, bamIndicies separate
    ifls <- c(tempfile(), tempfile())    # does not exist
    bv <- BamViews(fls, ifls)
    ifls <- setNames(ifls, names(bv))
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
    
.scanBam_ok <-
    function(target, current)
{
    checkIdentical(length(target), length(current))
    checkIdentical(class(target), class(current))
    for (i in seq_along(target)) {
        for (j in seq_along(target[[i]])) {
            t <- target[[i]][[j]]
            c <- current[[i]][[j]]
            switch(as.vector(class(t)),
                   PhredQuality=,
                   DNAStringSet={
                       checkIdentical(as.character(t), as.character(c))
                   },
                   checkIdentical(t, c))
        }
    }
}

test_BamViews_scanBam <- function()
{
    src <- system.file("unitTests", "cases", package="Rsamtools")
    fl <- c(system.file("extdata", "ex1.bam", package="Rsamtools"),
            file.path(src, "ex1_shuf1000.bam"))
    bv <- BamViews(fl)
    res <- scanBam(bv)
    for (i in seq_along(fl))
        .scanBam_ok(scanBam(fl[[i]]), res[[i]])

    param <- ScanBamParam(what="rname")
    res <- scanBam(bv, param=param)
    for (i in seq_along(fl))
        .scanBam_ok(scanBam(fl[[i]], param=param), res[[i]])
}

test_BamViews_countBam <- function()
{
    src <- system.file("unitTests", "cases", package="Rsamtools")
    fl <- c(system.file("extdata", "ex1.bam", package="Rsamtools"),
            file.path(src, "ex1_shuf1000.bam"))
    bv <- BamViews(fl)
    res <- countBam(bv)
    for (i in seq_along(fl))
        checkIdentical(countBam(fl[i]), res[[i]])

    which <- IRangesList(seq1=IRanges(1, 1000),
                         seq2=IRanges(1, 1000))
    param <- ScanBamParam(which=which)
    bamRanges <- GRanges(c("seq1", "seq2"), IRanges(c(1,1), 1000))
    bv <- BamViews(fl, bamRanges=bamRanges)
    res <- countBam(bv)
    for (i in seq_along(fl))
        checkIdentical(countBam(fl[i], param=param), res[[i]])
}
