setGeneric(".RsamtoolsFileList",
           function(file, ..., yieldSize=NA_integer_, class)
               standardGeneric(".RsamtoolsFileList"),
           signature="file")

setMethod(.RsamtoolsFileList, "missing",
    function(file, ..., classDef = class, yieldSize=NA_integer_, class)
{
    new(paste0(class, "List"))
})

setMethod(.RsamtoolsFileList, "character",
    function(file, index, ..., classDef=class, yieldSize=NA_integer_, class)
{
    fun <- function(elt, ..., yieldSize, classDef)
        do.call(classDef, list(elt, ..., yieldSize=yieldSize))
    if (is.null(names(file)))
        names(file) <- basename(file)
    listData <- if (!missing(index) && length(index))
        Map(fun, file, as.character(index), ...,
            MoreArgs=list(yieldSize=yieldSize, classDef=classDef))
    else if (missing(index))
        Map(fun, file, ...,
            MoreArgs=list(yieldSize=yieldSize, classDef=classDef)
        )
    else
        ## support old index=character() variant
        Map(fun, file, ..., MoreArgs=list(index=index,
                              yieldSize=yieldSize, classDef=classDef))
    new(paste0(class, "List"), listData=listData)
})

setMethod(.RsamtoolsFileList, "ANY",
    function(file, ..., classDef = class, yieldSize=NA_integer_, class)
{
    list <- list(file, ...)
    if (length(list) == 1 && (is.list(list[[1L]]) || is(list[[1L]], "List")))
        list <- as.list(list[[1L]])
    new(paste0(class, "List"), listData=list)
})

setMethod(.RsamtoolsFileList, "RsamtoolsFile",
    function(file, ..., classDef = class, yieldSize=NA_integer_, class)
{
    new(paste0(class, "List"), listData=list(file, ...))
})

setMethod(path, "RsamtoolsFileList",
    function(object, ...)
{
    vapply(object, path, character(1))
})

setMethod(index, "RsamtoolsFileList",
    function(object, ...)
{
    sapply(object, index, ...)
})

setReplaceMethod("index", "RsamtoolsFileList",
    function(object, ..., value)
{
    stopifnot(length(value) == length(path(object)))
    for (i in seq_along(object))
        index(object[[i]]) <- value[i]
    object
})

setMethod(yieldSize, "RsamtoolsFileList",
    function(object, ...)
{
    vapply(object, yieldSize, numeric(1))
})

setReplaceMethod("yieldSize", "RsamtoolsFileList", 
    function(object, ..., value)
{
    for (i in seq_along(object))
      yieldSize(object[[i]]) <- value
    object
})

setMethod(isOpen, "RsamtoolsFileList", 
    function(con, rw="") 
{
    sapply(as.list(con), isOpen, rw="read")
})

open.RsamtoolsFileList <-
    function(con, ...)
{
    for (f in as.list(con))
        open(f, ...)
    con
}

close.RsamtoolsFileList <-
    function(con, ...)
{
    for (f in as.list(con))
        close(f, ...)
    con
}

setMethod(names, "RsamtoolsFileList",
    function(x)
{
    nms <- callNextMethod()
    if (is.null(nms))
        nms <- sapply(x, function(elt) basename(path(elt)))
    nms
})

## implementations

BcfFileList <- function(...) .RsamtoolsFileList(..., class="BcfFile")

TabixFileList <- function(...)
    .RsamtoolsFileList(..., class="TabixFile")

FaFileList <- function(...) .RsamtoolsFileList(..., class="FaFile")

## BamFileList

setMethod(countBam, "BamFileList",
    function(file, index=file, ..., param=ScanBamParam())
{
    counts <- lapply(file, countBam, ..., param=param)
    do.call(rbind, counts)
})

setMethod(mergeBam, "BamFileList",
    function(files, destination, ...)
{
    files <- sapply(files, path)
    mergeBam(files, destination, ...)
})
