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
    cat("class:", class(object), "\n")
    cat(.ppath("path", path(object)))
    cat(.ppath("index", index(object)))
    cat("isOpen:", isOpen(object), "\n")
})
