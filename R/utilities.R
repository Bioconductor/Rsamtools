.ppath <- function(tag, filepath)
{
    wd <- options('width')[[1]] - nchar(tag) - 6
    if (0L == length(filepath) || nchar(filepath) < wd)
        return(sprintf("%s: %s\n", tag, filepath))
    bname <- basename(filepath)
    wd1 <- wd - nchar(bname)
    dname <- substr(dirname(filepath), 1, wd1)
    sprintf("%s: %s...%s%s\n",
            tag, dname, .Platform$file.sep, bname)
}

.io_check_exists <-
    function(file)
{
    idx <- !grepl("^(ftp)|(http)://", file)
    if (!all(sapply(file[idx], file.exists))) {
        msg <- paste(sprintf("'%s'", file[idx]), collapse="\n  ")
        stop("file(s) do not exist:\n  ", msg)
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
        file.rename(from, to) ||
            (file.copy(from, to) && file.remove(from))
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
    ok
}

.uunlist <-
    function(x) unlist(x, use.names=FALSE)

setMethod(.asSpace, "RangesList", function(x) {
    list(as.character(space(x)), .uunlist(start(x)), .uunlist(end(x)))
})

setMethod(.asSpace, "RangedData", function(x) .asSpace(ranges(x)))

setMethod(.asSpace, "GRanges", function(x) {
    list(as.character(seqnames(x)), start(x), end(x))
})

