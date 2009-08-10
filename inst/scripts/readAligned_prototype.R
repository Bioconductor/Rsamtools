library(Rsamtools)
library(ShortRead)
.file_names <- ShortRead:::.file_names

.readAligned_bam <-
    function(dirPath, pattern=character(0), ..., param=ScanBamParam())
{
    files <- .file_names(dirPath, pattern)
    param <-
        initialize(param,
                   what=c("qname", "flag", "rname", "strand", "pos",
                          "mapq", "cigar", "seq", "qual"))

    ## handle multiple files
    result <- lapply(files, scanBam, param=param, ...)
    if (length(bamWhich(param)) != 0)
        result <- unlist(result, recursive=FALSE, use.names=FALSE)

    ulist <- function(X, ...)
        unlist(lapply(X, "[[", ...), use.names=FALSE)
    uxscat <- function(X, ...)
        do.call(xscat, ulist(X, ...))

    ## FIXME: currently we don't represent these features:
    cigar <- ulist(result, "cigar")
    keep <- !grepl("[IDHSP]", cigar)
    if (!all(keep))
        warning("ignoring reads with any of I, D, H, S, P in 'cigar'")

    AlignedRead(sread=uxscat(result, "seq")[keep],
                id=BStringSet(ulist(result, "qname"))[keep],
                quality=FastqQuality(uxscat(result, "qual"))[keep],
                chromosome=factor(ulist(result, "rname"))[keep],
                strand=ulist(result, "strand")[keep],
                position=ulist(result, "pos")[keep],
                alignQuality=NumericQuality(ulist(result, "mapq"))[keep],
                alignData=AlignedDataFrame(
                  data=data.frame(
                    flag=ulist(result, "flag")[keep],
                    cigar=ulist(result, "cigar")[keep]),
                  metadata=data.frame(
                    labelDescription=c(
                      "Type of read; see ?scanBam",
                      "Alignement description, see ?scanBam"))))
}

##

fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
res <- .readAligned_bam(fl)
