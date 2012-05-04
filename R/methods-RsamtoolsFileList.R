setGeneric(".RsamtoolsFileList",
           function(..., yieldSize=NA_integer_, class)
               standardGeneric(".RsamtoolsFileList"),
           signature="...")

setMethod(.RsamtoolsFileList, "character",
    function(..., yieldSize=NA_integer_, class)
{
    listData <- lapply(..1, function(elt, yieldSize, class) {
        do.call(class, list(elt, yieldSize=yieldSize))
    }, yieldSize, class)
    new(paste(class, "List", sep=""), listData=listData)
})

setMethod(.RsamtoolsFileList, "ANY",
    function(..., yieldSize=NA_integer_, class)
{
    new(paste(class, "List", sep=""), listData=list(...))
})

setMethod(path, "RsamtoolsFileList",
    function(object, ...)
{
    sapply(as.list(object), path)
})

setMethod(yieldSize, "RsamtoolsFileList",
    function(object, ...)
{
    sapply(as.list(object), yieldSize)
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

BamFileList <-
    function(..., yieldSize=NA_integer_)
{
    .RsamtoolsFileList(..., yieldSize=yieldSize, class="BamFile")
}

BcfFileList <- function(...) .RsamtoolsFileList(..., class="BcfFile")

TabixFileList <-
    function(..., yieldSize=NA_integer_)
{
    .RsamtoolsFileList(..., yieldSize=yieldSize, class="TabixFile")
}

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
