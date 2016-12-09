.extptr <- function(object) object$.extptr

setMethod(index, "RsamtoolsFile",
    function(object, ..., asNA=TRUE)
{
    index <- object$index
    if (asNA && ((length(index) == 0L) || !nzchar(index)))
        NA_character_
    else
        index
})

setReplaceMethod("index", "RsamtoolsFile",
    function(object, ..., value)
{
    stopifnot(length(value) == 1L)
    object$index <- as.character(value)
    object
})

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

.RsamtoolsFile <-
    function(g, path, index, ..., yieldSize=NA_integer_)
{
    if (1L != length(path))
        stop("'file' must be length 1") # argh! public api is 'file'
    if (1L < length(index))
        stop("'index' must be length 0 or 1")
    if (1L != length(yieldSize))
        stop("'yieldSize' must be length 1")
    yieldSize <- as.integer(yieldSize)
    if (!(yieldSize > 0L || is.na(yieldSize)))
        stop("'yieldSize' must be >0 or NA")
    if (length(index) && is.na(index))
        index <- character(0)
    g$new(path=.normalizePath(path), index=.normalizePath(index), ...,
          yieldSize=yieldSize)
}

setMethod(path, "RsamtoolsFile", function(object, ...) object$path)

setMethod(isOpen, "RsamtoolsFile", function(con, rw="") FALSE)

setMethod(show, "RsamtoolsFile", function(object) {
    cat("class:", class(object), "\n")
    cat(.ppath("path", path(object)))
    cat(.ppath("index", index(object)))
    cat("isOpen:", isOpen(object), "\n")
    cat("yieldSize:", yieldSize(object), "\n")
})
