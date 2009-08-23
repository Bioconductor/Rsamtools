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
