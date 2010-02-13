library(Rsamtools)
library(ShortRead)
.file_names <- ShortRead:::.file_names
.throw <- ShortRead:::.throw

.readAligned_bamWhat <- function()
{
    c("qname", "flag", "rname", "strand", "pos", "mapq", "seq",
      "qual")
}

.readAligned_bam <-
    function(dirPath, pattern=character(0), ..., 
             param=ScanBamParam(
               simpleCigar=TRUE,
               reverseComplement=TRUE,
               what=.readAligned_bamWhat()))
{
    files <- 
        if (!grepl("^(ftp|http)://", dirPath))
            .file_names(dirPath, pattern)
        else {
            if (length(dirPath) != 1 || length(pattern) != 0) {
                msg <- paste("ftp:// and http:// support requires",
                             "'dirPath' as character(1),",
                             "'pattern' as character(0)", collapse="")
                .throw(SRError("UserArgumentMismatch", msg))
            }
            dirPath
        }
    ## FIXME: currently we only deal with cigars without indels
    if (!missing(param)) {
        if (bamSimpleCigar(param) != TRUE) {
            msg <- paste("using 'TRUE' for 'bamSimpleCigar(param)'",
                         "(skipping reads with I, D, H, S, or P in 'cigar')")
            msg <- paste(strwrap(msg, exdent=2), collapse="\n")
            warning(msg)
        }
        if (bamReverseComplement(param) != TRUE) {
            msg <- "using 'TRUE' for 'bamReverseComplement(param)'"
            msg <- paste(strwrap(msg, exdent=2), collapse="\n")
            warning(msg)
        }
        if (!setequal(bamWhat(param), .readAligned_bamWhat()))
        {
            msg <- sprintf("using '%s' for 'bamWhat(param)'",
                           paste(.readAligned_bamWhat(),
                                 collapse="', '"))
            msg <- paste(strwrap(msg, exdent=2), collapse="\n")
            warning(msg)
        }
        param <- initialize(param, simpleCigar=TRUE,
                            reverseComplement=TRUE,
                            what=.readAligned_bamWhat())
    }
    ## handle multiple files
    if (length(files) != 1)
        stop("dirPath, pattern must match exactly 1 file, but matched ",
             length(files))
##     result <- lapply(files, scanBam, param=param, ...)
    result <- scanBam(files, ..., param=param)
##     if (length(bamWhich(param)) != 0)
##         result <- unlist(result, recursive=FALSE, use.names=FALSE)

    ulist <- function(X, ...)
        unlist(lapply(X, "[[", ...), use.names=FALSE)
    uxsappend <- function(X, ...) {
        elts <- ulist(X, ...)
        res <- elts[[1]]
        for (i in seq_along(elts[-1]))
             res <- append(res, elts[[i+1]])
        res
    }

    AlignedRead(sread=uxsappend(result, "seq"),
                id=BStringSet(ulist(result, "qname")),
                quality=FastqQuality(as(uxsappend(result, "qual"),
                  "BStringSet")),
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

