.io_check_exists <-
    function(file)
{
    idx <- !grepl("^(ftp)|(http)://", file)
    if (!all(sapply(file[idx], file.exists))) {
        msg <- paste(sprintf("'%s'", file[idx]), collapse="\n  ")
        stop("'file' elements do not exist:\n  ", msg)
    }
}

.show_classname <-
    function(x) cat("class: ", class(x), "\n", sep="")

.normalizePath <-
    function(path)
{
    idx <- !grepl("^(ftp)|(http)://", path)
    path[idx] <- normalizePath(path.expand(path[idx]), mustWork=FALSE)
    path
}

.uunlist <-
    function(x) unlist(x, use.names=FALSE)
