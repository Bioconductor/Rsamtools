library(Rsamtools)

fl <- system.file("extdata", "ex1.bam", package="Rsamtools")

res <- scanBam(fl)
names(res)
lapply(res, head, 3)
table(res[["cigar"]])
table(width(res[["seq"]]))
table(res[["strand"]])
table(res[["flag"]])

which <- RangesList(seq1=IRanges(1000, 2000),
                    seq2=IRanges(c(100, 1000), c(1000, 2000)))
p1 <- ScanBamParam(which=which)
res <- scanBam(fl, param=p1)
names(res)
lapply(res[[1]], head, 3)

p2 <- ScanBamParam(what=c("rname", "strand", "pos", "width"))
lapply(scanBam(fl, param=p2), head, 3)
                
p3 <- ScanBamParam(flag=scanBamFlag(isPaired=TRUE),
                   what=c("rname", "strand", "pos", "width"))
lapply(scanBam(fl, param=p3), head, 3)

p4 <- ScanBamParam(flag=scanBamFlag(isMinusStrand=FALSE))
lapply(scanBam(fl, param=p4), head, 3)

.Call(Rsamtools:::.read_bam_header, fl, "rbu", TRUE)

## readAligned(fl, type="bam", which=which)

