setwd("/home/mtmorgan/proj/a/1000g/slx_maq_09_index_files")
ftpBase = "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/pilot_data/data/"
if (length(list.files(pattern=".*bai")) == 0)
{
    library(RCurl)
    indivs <- strsplit(getURL(ftpBase, ftplistonly=TRUE), "\n")[[1]]
    aln <- paste(ftpBase, indivs, "/alignment/", sep="")

    urls <- strsplit(getURL(aln, dirlistonly=TRUE), "\n")

    urls0 <- urls[sapply(urls, length) != 0]

    fls0 <- unlist(unname(urls0))
    fls1 <- fls0[grepl("bai$", fls0)]
    fls <- fls1[sapply(strsplit(fls, "\\."), length)==7]
    m <- t(as.data.frame(strsplit(fls, "\\.")))[,1:5]
    dimnames(m) <- list(NULL,
                        c("Individual", "Platform", "Alignment", "SRA", "Date"))
    df <- cbind(as.data.frame(m), File=fls)
    xtabs(~Platform+Alignment+Date, df)
    ## > xtabs(~Platform+Alignment, df)
    ##         Alignment
    ## Platform corona maq MOSAIK ssaha
    ##    454        0   0      9    14
    ##    SLX        0  25     11     0
    ##    SOLID      7   0      0     0

    urls1 <- Filter(function(x) length(x) != 0,
                    lapply(urls, function(x) x[grepl("SLX.maq.*2009_08.*bai$", x)]))
    slxMaq09 <- mapply(paste, names(urls1), urls1, sep="", USE.NAMES=FALSE)

    mapply(download.file, slxMaq09, basename(slxMaq09), MoreArgs=list(method="curl"))
}    
indexFiles <- list.files(pattern="bai$")
slxMaq09 <- paste(ftpBase, sub("^([^\\.]+).*$", "\\1", indexFiles),
                  "/alignment/", sub("\\.bai$", "", indexFiles), sep="")
## headers <- scanBamHeader(fls)           # nothing useful here :(

## Some regions of interest: genes involved in caffeine metabolism
library(KEGG.db)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg18)
library(biomaRt)
library(GenomicFeatures)
kid <- revmap(KEGGPATHID2NAME)[["Caffeine metabolism"]]
egid <- KEGGPATHID2EXTID[[sprintf("hsa%s", kid)]]

mart <- useMart("ensembl", "hsapiens_gene_ensembl")
ensid <- getBM(c("ensembl_transcript_id"), filters="entrezgene",
               values=egid, mart=mart)[[1]]
txdb <- makeTranscriptDbFromBiomart(transcript_ids=ensid)
saveFeatures(txdb, "caffeine-txdb.sqlite")

egid <- egid[egid != "9"]               # multiple locations
tbl <- merge(toTable(org.Hs.egCHRLOC[egid]), toTable(org.Hs.egCHRLOCEND[egid]))
rng <- with(tbl, {
    lvls <- sprintf("chr%s", sort(as.integer(unique(Chromosome))))
    rng0 <- GRanges(Chromosome,
                    IRanges(start=abs(start_location), end=abs(end_location)),
                    strand=ifelse(start_location>=0, "+", "-"),
                    EntrezId=gene_id,
                    Genename=mappedRkeys(org.Hs.egGENENAME[gene_id]))
    metadata(rng0) <- list(Genome="hg18?")
    rng0
})

library(GenomicAlignments)
library(multicore)
bv <- BamViews(fls, sub(".bai", "", indexFiles),
               bamRanges=GRanges(seqnames=seqnames(rng),
                 ranges=ranges(rng),
                 EntrezId=mcols(rng)[["EntrezId"]]))
gapped <- readGAlignmentsFromBam(bv[1:2,1:3])
