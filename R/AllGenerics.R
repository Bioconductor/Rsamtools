setGeneric("scanBam",
           function(file, index=file, ...) standardGeneric("scanBam"))

setGeneric("scanBamHeader",
           function(files, ...) standardGeneric("scanBamHeader"))

setGeneric("countBam",
           function(file, index=file, ...) standardGeneric("countBam"))

setGeneric("filterBam",
           function(file, destination, index=file, ...)
           standardGeneric("filterBam"))

setGeneric("indexBam",
           function(files, ...) standardGeneric("indexBam"))

setGeneric("readPileup",
           function(file, ...) standardGeneric("readPileup"))

setGeneric("readBamGappedAlignments",
           function(file, index, ..., which, ans.subtype="Alignments0")
               standardGeneric("readBamGappedAlignments"),
           signature="file")

