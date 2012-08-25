.ScanBcfParam <-
    function(fixed=character(), info=character(), geno=character(), 
             trimEmpty=TRUE, which, class="ScanBcfParam")
{
    if (1L == length(fixed) && is.na(fixed))
        fixed <- as.character(fixed)
    if (1L == length(info) && is.na(info))
        info <- as.character(info)
    if (1L == length(geno) && is.na(geno))
        geno <- as.character(geno)
    new(class, which=which, fixed=fixed, info=info, geno=geno,
        trimEmpty=trimEmpty)
}

## ScanBcfParam

setMethod(ScanBcfParam, c(which="missing"),
    function(fixed=character(), info=character(), geno=character(), 
             trimEmpty=TRUE, which, ...)
{
    which <- IRangesList()
    names(which) <- character()
    .ScanBcfParam(fixed, info, geno, trimEmpty, which, ...)
})

setMethod(ScanBcfParam, c(which="RangesList"), 
    function(fixed=character(), info=character(), geno=character(), 
             trimEmpty=TRUE, which, ...)
{
    .ScanBcfParam(fixed, info, geno, trimEmpty, which, ...)
})

setMethod(ScanBcfParam, c(which="RangedData"),
    function(fixed=character(), info=character(), geno=character(), 
             trimEmpty=TRUE, which, ...)
{
    which <- ranges(which)
    .ScanBcfParam(fixed, info, geno, trimEmpty, which, ...)
})

setMethod(ScanBcfParam, c(which="GRanges"),
    function(fixed=character(), info=character(), geno=character(), 
             trimEmpty=TRUE, which, ...)
{
    which <- split(ranges(which), seqnames(which))
    .ScanBcfParam(fixed, info, geno, trimEmpty, which, ...)
})

## accessors

bcfFixed <- function(object) slot(object, "fixed")
bcfInfo <- function(object) slot(object, "info")
bcfGeno <- function(object) slot(object, "geno")
bcfTrimEmpty <- function(object) slot(object, "trimEmpty")
bcfWhich <- function(object) slot(object, "which")

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
        paste0(cl, lbl)
    }
    cat("class:", class(object), "\n")
    cat(sprintf("%s: %d elements\n", .clslbl("Which"),
                length(bcfWhich(object))))
    cat(.clslbl("Fixed:"), .ptags(bcfFixed(object)), "\n")
    cat(.clslbl("Info:"), .ptags(bcfInfo(object)), "\n")
    cat(.clslbl("Geno:"), .ptags(bcfGeno(object)), "\n")
})

