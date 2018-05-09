bgzipTabix <-
    function(fromFname, toFname = paste(fromFname, "gz", sep="."),
             overwrite=FALSE)
{
    .Defunct("bgzip", package="Rsamtools")
}

.zip <-
    function(func, file, dest, overwrite)
{
    file <- .normalizePath(file)
    dest <- .normalizePath(dest)
    if (!is.character(dest) || 1L != length(dest))
        stop("'dest' must be character(1)")
    if (!overwrite && file.exists(dest))
        stop("'dest' exists:\n  dest: ", dest)
    tryCatch({
        .Call(func, file, dest)
    }, error=function(err) {
        msg <- sprintf("'%s' error: %s\n  file: %s\n  dest: %s",
                       sub(".", "", func, fixed=TRUE), conditionMessage(err),
                       file, dest)
        stop(msg, call.=FALSE)
    })
}

bgzip <-
    function(file, dest = sprintf("%s.bgz", sub("\\.gz$", "", file)),
             overwrite=FALSE)
{
    .zip(.bgzip, file, dest, overwrite)
}


razip <-
    function(file, dest = sprintf("%s.rz", sub("\\.gz$", "", file)),
             overwrite=FALSE)
{
    .Defunct("bgzip")
}
