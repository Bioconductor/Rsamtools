.readPileup_table <-
    function(conn, colClasses, ...)
{
    read.table(conn, colClasses=colClasses,
               col.names=names(colClasses), sep="\t", header=FALSE,
               quote="", comment.char="", fill=TRUE, ...)
}

.readPileup_indel_idx <- function(df) which(df[[3]] == "*")

.readPileup_SNP <-
    function(file, ..., variant)
{
    colClasses <-
        c(space="factor", position="integer",
          referenceBase="character", consensusBase="character",
          consensusQuality="integer", snpQuality="integer",
          maxMappingQuality="integer", coverage="integer", "NULL",
          "NULL", "NULL", "NULL", "NULL", "NULL", "NULL")
    dat <- .readPileup_table(file, colClasses, ...)
    idx <- .readPileup_indel_idx(dat)
    if (length(idx) > 0L) {
        if (variant == "SNP")
            idx <- c(idx, idx-1)
        dat <- dat[-idx,]
    }
    GRanges(seqnames=dat[,1],
            ranges=IRanges(start=dat[,2],end=dat[,2]),
            referenceBase=factor(dat[,3], levels=DNA_ALPHABET),
            consensusBase=factor(dat[,4], levels=DNA_ALPHABET),
            consensusQuality=dat[,5],
            snpQuality=dat[,6],
            maxMappingQuality=dat[,7],
            coverage=dat[,8])
}

.readPileup_indel <-
    function(file, ...)
{
    colClasses <- c(space="factor", position="integer",
                    reference="character", consensus="character",
                    consensusQuality="integer", snpQuality="integer",
                    maxMappingQuality="integer", coverage="integer",
                    alleleOne="character", alleleTwo="character",
                    alleleOneSupport="integer", alleleTwoSupport="integer",
                    additionalIndels="integer", "NULL", "NULL")
    dat <- .readPileup_table(file, colClasses, ...)
    idx <- .readPileup_indel_idx(dat)
    if (length(idx) != 0L) {
        dat0 <- dat[idx-1,]
        dat <- dat[idx,]
    } else {
        dat <- dat0 <- dat[FALSE,]
    }

    GRanges(seqnames=dat[,1],
            ranges=IRanges(start=dat0[,2],end=dat0[,2]),
            referenceBase=factor(dat0[,3], levels=DNA_ALPHABET),
            consensusBase=factor(dat0[,4]),
            consensusQuality=dat0[,5],
            snpQuality=dat0[,6],
            maxMappingQuality=dat0[,7],
            coverage=dat0[,8],
            alleleOne=dat[,9],
            alleleOneSupport=dat[,11],
            alleleTwo=dat[,10],
            alleleTwoSupport=dat[,12],
            additionalIndels=dat[,13])
}

setMethod(readPileup, "connection",
          function(file, ..., variant=c("SNP", "indel", "all"))
{
    variant <- match.arg(variant)
    switch(variant, SNP=, all=.readPileup_SNP(file=file, ...,
           variant=variant), indel=.readPileup_indel(file=file, ...))
})

setMethod(readPileup, "character", function(file, ...)
{
    conn <- file(file, "r")
    on.exit(close(conn))
    callGeneric(conn, ...)
})
