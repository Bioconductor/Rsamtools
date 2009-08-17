library(Rsamtools)
library(ShortRead)
.file_names <- ShortRead:::.file_names

.readAligned_bamWhat <- function()
{
    c("qname", "flag", "rname", "strand", "pos", "mapq", "seq",
      "qual")
}

.readAligned_bam <-
    function(dirPath, pattern=character(0), ..., 
             param=ScanBamParam(
               simpleCigar=TRUE,
               what=.readAligned_bamWhat()))
{
    files <- .file_names(dirPath, pattern)
    ## FIXME: currently we only deal with cigars without indels
    if (!missing(param)) {
        if (bamSimpleCigar(param) != TRUE) {
            msg <- paste("using 'TRUE' for 'bamSimpleCigar(param)'",
                         "(skipping reads with I, D, H, S, or P in 'cigar')")
            msg <- paste(strwrap(msg, exdent=2),
                         collapse="\n")
            warning(msg)
        }
        if (!setequal(bamWhat(param), .readAligned_bamWhat()))
        {
            msg <- sprintf("using '%s' for 'bamWhat(param)'",
                           paste(.readAligned_bamWhat(),
                                 collapse="', '"))
            msg <- paste(strwrap(msg, exdent=2),
                         collapse="\n")
            warning(msg)
        }
        param <- initialize(param, simpleCigar=TRUE,
                            what=.readAligned_bamWhat())
    }
    ## handle multiple files
    result <- lapply(files, scanBam, param=param, ...)
    if (length(bamWhich(param)) != 0)
        result <- unlist(result, recursive=FALSE, use.names=FALSE)

    ulist <- function(X, ...)
        unlist(lapply(X, "[[", ...), use.names=FALSE)
    uxscat <- function(X, ...)
        do.call(xscat, ulist(X, ...))

    AlignedRead(sread=uxscat(result, "seq"),
                id=BStringSet(ulist(result, "qname")),
                quality=FastqQuality(uxscat(result, "qual")),
                chromosome=factor(ulist(result, "rname")),
                strand=ulist(result, "strand"),
                position=ulist(result, "pos"),
                alignQuality=NumericQuality(ulist(result, "mapq")),
                alignData=AlignedDataFrame(
                  data=data.frame(
                    flag=ulist(result, "flag")),
                  metadata=data.frame(
                    labelDescription=c(
                      "Type of read; see ?scanBam"))))
}

##

fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
res <- .readAligned_bam(fl, param=ScanBamParam())
res <- .readAligned_bam(fl, param=ScanBamParam(simpleCigar=TRUE))
(res <- .readAligned_bam(fl))
