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

setGeneric("rname", function(x) standardGeneric("rname"))

setGeneric("rname<-", function(x, value) standardGeneric("rname<-"))

setGeneric("cigar", function(x) standardGeneric("cigar"))

setGeneric("qwidth", function(x) standardGeneric("qwidth"))

setGeneric("grglist", function(x) standardGeneric("grglist"))

setGeneric("grg", function(x) standardGeneric("grg"))

setGeneric("rglist", function(x) standardGeneric("rglist"))

setGeneric("readBAMasAlignments0",
           function(file, index, ..., which)
               standardGeneric("readBAMasAlignments0"),
           signature="file")

setGeneric("readBAMasAlignments1",
           function(file, index, ..., which)
               standardGeneric("readBAMasAlignments1"),
           signature="file")

setGeneric("readBAMasAlignments2",
           function(file, index, ..., which)
               standardGeneric("readBAMasAlignments2"),
           signature="file")

setGeneric("updateCigarAndStart",
           function(x, cigar=NULL, start=NULL)
               standardGeneric("updateCigarAndStart"))

setGeneric("qnarrow",
           function(x, start=NA, end=NA, width=NA)
               standardGeneric("qnarrow"))

