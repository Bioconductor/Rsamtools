.ScanBcfParam <-
    function(info=character(), geno=character(), trimEmpty=TRUE, which)
{
    if (1L == length(info) && is.na(info))
        info <- as.character(info)
    if (1L == length(geno) && is.na(geno))
        geno <- as.character(geno)
    new("ScanBcfParam", which=which, info=info, geno=geno,
        trimEmpty=trimEmpty)
}

setMethod(ScanBcfParam, c(which="missing"),
    function(info=character(), geno=character(), trimEmpty=TRUE, which)
{
    which <- IRangesList()
    names(which) <- character()
    .ScanBcfParam(info, geno, trimEmpty, which)
})

setMethod(ScanBcfParam, c(which="RangesList"), .ScanBcfParam)

setMethod(ScanBcfParam, c(which="RangedData"),
    function(info=character(), geno=character(), trimEmpty=TRUE, which)
{
    which <- ranges(which)
    .ScanBcfParam(info, geno, trimEmpty, which)
})

setMethod(ScanBcfParam, c(which="GRanges"),
    function(info=character(), geno=character(), trimEmpty=TRUE, which)
{
    which <- split(ranges(which), seqnames(which))
    .ScanBcfParam(info, geno, trimEmpty, which)
})

bcfInfo <- function(object) slot(object, "info")
bcfGeno <- function(object) slot(object, "geno")
bcfTrimEmpty <- function(object) slot(object, "trimEmpty")
bcfWhich <- function(object) slot(object, "which")

setMethod(show, "ScanBcfParam", function(object) 
{
    .ptags <- function(tags) {
        if (length(tags))
            paste(tags, collapse=", ")
        else "character() [All]"
    }
    cat("class;", class(object), "\n")
    cat(sprintf("bcfWhich: %d elements\n", length(bcfWhich(object))))
    cat("bcfInfo:", .ptags(bcfInfo(object)), "\n")
    cat("bcfGeno:", .ptags(bcfGeno(object)), "\n")
})
