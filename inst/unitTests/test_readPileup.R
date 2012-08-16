pileupFile <- system.file("extdata", "pileup.txt", package="Rsamtools")
src <- system.file("unitTests", "cases", package="Rsamtools")
nostarsFile <- file.path(src, "pileup-no-stars.txt")


.check_SNP_pileup <- function(tbl, rdf)
{
    checkIdentical(nrow(tbl), length(rdf))
    checkIdentical(as.character(tbl[[1]]), as.character(seqnames(rdf)))
    checkIdentical(tbl[[2]], start(rdf))
    for (i in 3:8)
        checkIdentical(as.character(tbl[[i]]),
                       as.character(mcols(rdf)[[i-2]]))
}

test_readPileup <- function()
{
    tbl <- read.table(pileupFile, fill=TRUE, quote="", comment="",
                      col.names=1:16)

    rdf <- readPileup(pileupFile, variant="indel")
    idx <- which(tbl[[3]]=="*")
    checkIdentical(length(idx), length(rdf))

    rdf <- readPileup(pileupFile, variant="SNP")
    checkIdentical(nrow(tbl) - 2L * length(idx), length(rdf))
    .check_SNP_pileup(tbl[-c(idx, idx-1),], rdf)
}

test_readPileup_nostars <- function()
{
    rdf <- readPileup(nostarsFile)
    tbl <- read.table(nostarsFile)
    .check_SNP_pileup(tbl, rdf)

    rdf <- readPileup(nostarsFile, variant="indel")
    checkIdentical(0L, length(rdf))
}
