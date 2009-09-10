fl0 = "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/pilot_data/data/NA19240/alignment/NA19240.chrom6.SLX.maq.SRP000032.2009_07.bam"

dir <- getwd()
fl <- file.path(dir, "NA19240.chrom6.SLX.maq.SRP000032.2009_07.bam")
stopifnot(file.exists(paste(fl, "bai", sep="."))) # require local index

library(GenomicFeatures)
library(IRanges)
data(geneHuman)
transcripts <- transcripts(geneHuman, proximal=300)

library(Rsamtools)
## remote access
chr6a0 <- ranges(transcripts)[["chr6"]][1:2]
p10 <- ScanBamParam(which=RangesList(`6`=chr6a0))
(cnt0 <- countBam(fl0, fl0, param=p10))
sum(cnt0$records)
system.time(res0 <- scanBam(fl0, fl0, param=p10))


## local
if (file.exists(fl)) {
    cnt <- countBam(fl, param=p10); sum(cnt$records)
    system.time(res0 <- scanBam(fl, param=p10))

    ## larger, local
    chr6a <- ranges(transcripts)[["chr6"]][1:50]
    p1 <- ScanBamParam(which=RangesList(`6`=chr6a))
    sum(countBam(fl, param=p1)$records)
    system.time(res <- scanBam(fl, param=p1)) # about 30s
}
