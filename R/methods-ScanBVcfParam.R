.ScanBcfParam <-
    function(info=character(), geno=character(), trimEmpty=TRUE, which,
             class="ScanBcfParam")
{
    if (1L == length(info) && is.na(info))
        info <- as.character(info)
    if (1L == length(geno) && is.na(geno))
        geno <- as.character(geno)
    new(class, which=which, info=info, geno=geno,
        trimEmpty=trimEmpty)
}

## ScanBcfParam

setMethod(ScanBcfParam, c(which="missing"),
    function(info=character(), geno=character(), trimEmpty=TRUE, which,
             ...)
{
    which <- IRangesList()
    names(which) <- character()
    .ScanBcfParam(info, geno, trimEmpty, which, ...)
})

setMethod(ScanBcfParam, c(which="RangesList"), 
    function(info=character(), geno=character(), trimEmpty=TRUE, which,
             ...)
{
    .ScanBcfParam(info, geno, trimEmpty, which, ...)
})

setMethod(ScanBcfParam, c(which="RangedData"),
    function(info=character(), geno=character(), trimEmpty=TRUE, which,
             ...)
{
    which <- ranges(which)
    .ScanBcfParam(info, geno, trimEmpty, which, ...)
})

setMethod(ScanBcfParam, c(which="GRanges"),
    function(info=character(), geno=character(), trimEmpty=TRUE, which,
             ...)
{
    which <- split(ranges(which), seqnames(which))
    .ScanBcfParam(info, geno, trimEmpty, which, ...)
})

## ScanVcfParam

setMethod(ScanVcfParam, "ANY",
    function(info=character(), geno=character(), trimEmpty=TRUE, which, ...)
{
    ScanBcfParam(info, geno, trimEmpty, which, class="ScanVcfParam")
})

setMethod(ScanVcfParam, "missing",
    function(info=character(), geno=character(), trimEmpty=TRUE, which, ...)
{
    ScanBcfParam(info, geno, trimEmpty, class="ScanVcfParam")
})

## accessors

vcfInfo <- bcfInfo <- function(object) slot(object, "info")
vcfGeno <- bcfGeno <- function(object) slot(object, "geno")
vcfTrimEmpty <- bcfTrimEmpty <- function(object) slot(object, "trimEmpty")
vcfWhich <- bcfWhich <- function(object) slot(object, "which")

setMethod(show, "ScanBVcfParam", function(object) 
{
    .ptags <- function(tags) {
        if (length(tags))
            paste(tags, collapse=", ")
        else "character() [All]"
    }
    .clslbl <- function(lbl) {
        cl <- 
            if ("ScanBcfParam" == class(object)) "bcf" else "vcf"
        paste(cl, lbl, sep="")
    }
    cat("class:", class(object), "\n")
    cat(sprintf("%s: %d elements\n", .clslbl("Which"),
                length(bcfWhich(object))))
    cat(.clslbl("Info:"), .ptags(bcfInfo(object)), "\n")
    cat(.clslbl("Geno:"), .ptags(bcfGeno(object)), "\n")
})
