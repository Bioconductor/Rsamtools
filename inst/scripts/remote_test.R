suppressMessages(library(Rsamtools))
fl = "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/pilot_data/data/NA19240/alignment/NA19240.chrom6.SLX.maq.SRP000032.2009_07.bam"
p1 <- ScanBamParam(which=RangesList("6"=IRanges(10000, 11000)))
res <- scanBam(fl, param=p1)[[1]]
res[["seq"]]

fl = "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/pilot_data/data/NA19240/alignment/NA19240.chrom6.454.ssaha.SRP000032.2009_07.bam"
p1 <- ScanBamParam(which=RangesList("6"=IRanges(100000, 110000)))
res <- scanBam(fl, param=p1)[[1]]
res[["seq"]]

header <- scanBamHeader(fl)
txt <- header[[1]][[2]]
table(names(txt))
txt[names(txt) == "@HD"]                # 'header'
txt[names(txt) == "@RG"][1:5]           # 'read group'
txt[names(txt) == "@SQ"][1:5]           # 'sequence group'

