setMethod(scanBamHeader, "character", 
          function(files, ...)
{
    files <- .normalizePath(files)
    lst <- lapply(files, function(file) {
        bam <- openBam(file, character(0))
        on.exit(close(bam))
        header <- scanBamHeader(bam, ...)
        text <- strsplit(header[["text"]], "\n")[[1]]
        tag <- sub("^(@[A-Z]{2}).*", "\\1", text)
        text <- strsplit(sub("^@[A-Z]{2}\t(.*)", "\\1", text), "\t")
        names(text) <- tag
        header[["text"]] <- text
        header
    })
    names(lst) <- files
    lst
})
