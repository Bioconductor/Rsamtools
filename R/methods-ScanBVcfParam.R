.ScanBcfParam <-
    function(info=character(), geno=character(), trimEmpty=TRUE, which,
             asGRanges=character(), class="ScanBcfParam")
{
    if (1L == length(info) && is.na(info))
        info <- as.character(info)
    if (1L == length(geno) && is.na(geno))
        geno <- as.character(geno)
    if (1L == length(asGRanges) && is.na(asGRanges))
        asGRanges <- as.character(asGRanges)
    else
        stopifnot(asGRanges %in% c('info', 'geno'))

    new(class, which=which, info=info, geno=geno,
        trimEmpty=trimEmpty, asGRanges=asGRanges)
}

## ScanBcfParam

setMethod(ScanBcfParam, c(which="missing"),
    function(info=character(), geno=character(), trimEmpty=TRUE, which,
             asGRanges=character(), ...)
{
    which <- IRangesList()
    names(which) <- character()
    .ScanBcfParam(info, geno, trimEmpty, which, asGRanges, ...)
})

setMethod(ScanBcfParam, c(which="RangesList"), 
    function(info=character(), geno=character(), trimEmpty=TRUE, which,
             asGRanges=character(), ...)
{
    .ScanBcfParam(info, geno, trimEmpty, which, asGRanges, ...)
})

setMethod(ScanBcfParam, c(which="RangedData"),
    function(info=character(), geno=character(), trimEmpty=TRUE, which,
             asGRanges, ...)
{
    which <- ranges(which)
    .ScanBcfParam(info, geno, trimEmpty, which, asGRanges, ...)
})

setMethod(ScanBcfParam, c(which="GRanges"),
    function(info=character(), geno=character(), trimEmpty=TRUE, which,
             asGRanges, ...)
{
    which <- split(ranges(which), seqnames(which))
    .ScanBcfParam(info, geno, trimEmpty, which, asGRanges, ...)
})

## accessors

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
        paste(cl, lbl, sep="")
    }
    cat("class:", class(object), "\n")
    cat(sprintf("%s: %d elements\n", .clslbl("Which"),
                length(bcfWhich(object))))
    cat(.clslbl("Info:"), .ptags(bcfInfo(object)), "\n")
    cat(.clslbl("Geno:"), .ptags(bcfGeno(object)), "\n")
    cat(.clslbl("AsGRanges:"), .ptags(vcfAsGRanges(object)), "\n")
})

