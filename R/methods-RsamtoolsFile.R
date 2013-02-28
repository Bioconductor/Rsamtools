.extptr <- function(object) object$.extptr

index <- function(object) object$index

setMethod(yieldSize, "RsamtoolsFile",
    function(object, ...)
{
    object$yieldSize
})

setReplaceMethod("yieldSize", "RsamtoolsFile", 
    function(object, ..., value)
{
    if (1L != length(value))
        stop("'value' must be length 1")
    object$yieldSize <- as.integer(value)
    object
})

setMethod(obeyQname, "RsamtoolsFile",
    function(object, ...)
{
    object$obeyQname
})

setReplaceMethod("obeyQname", "RsamtoolsFile", 
    function(object, ..., value)
{
    if (1L != length(value))
        stop("'value' must be length 1")
    if (!is.logical(value))
        stop("'value' must be logical")
    object$obeyQname <- value
    object
})

.RsamtoolsFile <-
    function(g, path, index, ..., yieldSize=NA_integer_,
             obeyQname=FALSE)
{
    if (1L != length(yieldSize))
        stop("'yieldSize' must be length 1")
    yieldSize <- as.integer(yieldSize)
    if (!(yieldSize > 0L || is.na(yieldSize)))
        stop("'yieldSize' must be >0 or NA")
    g$new(path=.normalizePath(path), index=.normalizePath(index), ...,
          yieldSize=yieldSize, obeyQname=obeyQname, ...)
}

setMethod(path, "RsamtoolsFile", function(object, ...) object$path)

setMethod(isOpen, "RsamtoolsFile", function(con, rw="") FALSE)

setMethod(show, "RsamtoolsFile", function(object) {
    cat("class:", class(object), "\n")
    cat(.ppath("path", path(object)))
    cat(.ppath("index", index(object)))
    cat("isOpen:", isOpen(object), "\n")
    cat("yieldSize:", yieldSize(object), "\n")
    cat("obeyQname:", obeyQname(object), "\n")
})
