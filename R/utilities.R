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

.file.rename <-
    function(from, to)
{
    warn <- err <- NULL
    ok <- withCallingHandlers(tryCatch({
        file.rename(from, to)
    }, error=function(e) {
        err <<- append(err, conditionMessage(e))
        NULL
    }), warning=function(w) {
        warn <<- append(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
    }) || withCallingHandlers(tryCatch({
        file.copy(from, to) && file.remove(from)
    }, error=function(e) {
        err <<- append(err, conditionMessage(e))
        NULL
    }), warning=function(w) {
        warn <<- append(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    if (!ok) {
        msg <- "file.rename or file.copy/file.remove failed:\n  from: %s\n  to: %s\n  message(s): %s"
        stop(sprintf(msg, from, to, paste(c(warn, err), collapse="\n      ")))
    }
}

.uunlist <-
    function(x) unlist(x, use.names=FALSE)
