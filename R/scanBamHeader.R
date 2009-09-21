setMethod(scanBamHeader, "character", 
          function(files, ...)
{
    if (!is.character(files))
        stop("'file' must be character()")
    files <- .normalizePath(files)
    lst <- lapply(files, function(file) {
        header <- .Call(.read_bam_header, file, "rb")
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
