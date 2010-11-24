.extptr <- function(object) object$.extptr

path <- function(object) object$path

index <- function(object) object$index

.RsamtoolsFile <-
    function(g, path, index, ...)
{
    g$new(path=.normalizePath(path), index=.normalizePath(index), ...)
}

setMethod(isOpen, "RsamtoolsFile", function(con, rw="") FALSE)

setMethod(show, "RsamtoolsFile", function(object) {
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
    cat("class:", class(object), "\n")
    cat(.ppath("path", path(object)))
    cat(.ppath("index", index(object)))
    cat("isOpen:", isOpen(object), "\n")
})
