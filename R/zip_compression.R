bgzipTabix <-
    function(fromFname, toFname = paste(fromFname, "gz", sep="."),
             overwrite=FALSE)
{
    .Deprecated("bgzip", package="Rsamtools")
    bgzip(fromFname, toFname, overwrite)
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
    .Call(func, file, dest)
}

bgzip <-
    function(file, dest = sprintf("%s.gz", file), overwrite=FALSE)
{
    .zip(.bgzip, file, dest, overwrite)
}


razip <-
    function(file, dest = sprintf("%s.rz", file), overwrite=FALSE)
{
    .zip(.razip, file, dest, overwrite)
}
