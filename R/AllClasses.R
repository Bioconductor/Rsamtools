setOldClass("connection")
setOldClass(c("file", "connection"))
setOldClass(c("url", "connection"))
setOldClass(c("gzfile", "connection"))
setOldClass(c("bzfile", "connection"))
setOldClass(c("unz", "connection"))
setOldClass(c("file", "connection"))
setOldClass(c("pipe", "connection"))
setOldClass(c("fifo", "connection"))

setGeneric(".validity", function(object) standardGeneric(".validity"))

setClass("ScanBamParam",
         representation=representation(
           flag="integer",
           simpleCigar="logical",
           which="RangesList",
           what="character"))

setClass("Cigar",
         representation=representation("Sequence",
           .cigar="factor"),
         validity=.validity)
