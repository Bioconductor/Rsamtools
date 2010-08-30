source(system.file("scripts", "readAligned_prototype.R",
                   package="Rsamtools"))
fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
.what <- c('qname', 'flag', 'rname', 'strand', 'pos', 'mapq', 'seq',
           'qual')

test_readAligned_bam <- function()
{
    aln <- .readAligned_bam(fl)
    checkEquals(structure("AlignedRead", package = "ShortRead"),
                class(aln))
    checkTrue(validObject(aln))
    checkIdentical(3278L, length(aln))

    exp <- structure(c(6L, 37L, 2839L, 283L, 113L), .Dim = 5L,
                     .Dimnames = structure(list( c("33", "34", "35",
                     "36", "40")), .Names = ""), class = "table")
    checkIdentical(exp, table(width(aln)))

    exp <- structure(c(1498L, 1780L), .Dim = 2L, .Dimnames =
                     structure(list( c("seq1", "seq2")), .Names = ""),
                     class = "table")
    checkIdentical(exp, table(chromosome(aln)))

    exp <- structure(c(1631L, 1611L, 0L, 36L), .Dim = 4L, .Dimnames =
                     structure(list( c(levels(strand()), NA)), .Names =
                     ""), class = "table")
    checkIdentical(exp, table(strand(aln), useNA="always"))

    checkEquals(793.298, mean(position(aln), na.rm=TRUE),
                tolerance=10e-4)

    exp <- structure(c(3278L, 1L),
                     .Names = c("readName", "alignColumn"))
    checkIdentical(exp, dim(alignData(aln)))

    exp <- structure(c(35946L, 21903L, 21933L, 35608L, 139L), .Names =
                     c("A", "C", "G", "T", "other"))
    checkIdentical(exp,
                   alphabetFrequency(sread(aln), collapse=TRUE,
                                     baseOnly=TRUE))

    exp <- structure(c(263, 1, 20, 178, 575, 163, 195, 286, 285, 850,
                       290, 367, 416, 338, 391, 604, 600, 688, 894,
                       778, 1210, 2284, 1804, 1879, 3307, 10557,
                       79755, 6380, 150, 18, 3), .Names = c("!", "#",
                       "$", "%", "&", "'", "(", ")", "*", "+", ",",
                       "-", ".", "/", "0", "1", "2", "3", "4", "5",
                       "6", "7", "8", "9", ":", ";", "<", "=", ">",
                       "?", "@"))
    checkIdentical(exp, rowSums(consensusMatrix(quality(quality(aln)))))
}

.checkEquals0 <- function(aln0, aln1)
{
    checkIdentical(class(aln0), class(aln1))
    checkIdentical(length(aln0), length(aln1))
    checkIdentical(width(aln0), width(aln1))
    checkIdentical(chromosome(aln0), chromosome(aln1))
    checkIdentical(position(aln0), position(aln1))
    checkIdentical(strand(aln0), strand(aln1))
    checkIdentical(alignQuality(aln0), alignQuality(aln1))
    checkIdentical(as.character(sread(aln0)),
                   as.character(sread(aln1)))
    checkIdentical(as.character(quality(quality(aln0))),
                   as.character(quality(quality(aln1))))
}

test_readAligned_bam_isSimpleCigar <- function()
{
    p <- ScanBamParam(simpleCigar=TRUE, reverseComplement=TRUE,
                      what=.what)
    .checkEquals0(.readAligned_bam(fl, param=p),
                  .readAligned_bam(fl))
}

test_readAligned_bam_which <- function()
{
    which <- RangesList(seq1=IRanges(1000, 2000),
                        seq2=IRanges(c(100, 1000), c(1000, 2000)))
    p1 <- ScanBamParam(simpleCigar=TRUE, reverseComplement=TRUE,
                       what=.what, which=which)
    aln <- .readAligned_bam(fl, param=p1)
    checkEquals(structure("AlignedRead", package = "ShortRead"),
                class(aln))
    checkTrue(validObject(aln))

    res <- scanBam(fl, param=p1)
    exp <- sum(sapply(res, function(elt) length(elt[["pos"]])))
    checkEquals(exp, length(aln))

    exp <- structure(c(4L, 24L, 2076L, 205L, 88L), .Dim = 5L,
                     .Dimnames = structure(list( c("33", "34", "35",
                     "36", "40")), .Names = ""), class = "table")
    checkIdentical(exp, table(width(aln)))
}

