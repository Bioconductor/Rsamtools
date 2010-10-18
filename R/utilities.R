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
