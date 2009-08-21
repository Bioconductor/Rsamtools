source(system.file("scripts", "readAligned_prototype.R",
                   package="Rsamtools"))
fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
res <- .readAligned_bam(fl, param=ScanBamParam())
res <- .readAligned_bam(fl, param=ScanBamParam(simpleCigar=TRUE))
(res <- .readAligned_bam(fl))
sread(res)

which <- RangesList(seq1=IRanges(1000, 2000),
                    seq2=IRanges(c(100, 1000), c(1000, 2000)))
p1 <- ScanBamParam(which=which)
sread(.readAligned_bam(fl, param=p1))
