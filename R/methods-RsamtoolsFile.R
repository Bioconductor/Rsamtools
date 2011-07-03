.extptr <- function(object) object$.extptr

path <- function(object) object$path

index <- function(object) object$index

.RsamtoolsFile <-
    function(g, name, path, index, ...)
{
    if (1L != length(name))
        stop("'name' must be character(1)")
    g$new(name=name, path=.normalizePath(path),
          index=.normalizePath(index), ...)
}

setMethod(isOpen, "RsamtoolsFile", function(con, rw="") FALSE)

setMethod(names, "RsamtoolsFile", function(x) x$name)

setMethod(show, "RsamtoolsFile", function(object) {
    cat("class:", class(object), "\n")
    cat("names:", names(object), "\n")
    cat(.ppath("path", path(object)))
    cat(.ppath("index", index(object)))
    cat("isOpen:", isOpen(object), "\n")
})
