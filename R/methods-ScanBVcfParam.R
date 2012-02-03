.ScanBcfParam <-
    function(info=character(), geno=character(), trimEmpty=TRUE, which,
             asGRanges=FALSE, class="ScanBcfParam")
{
    if (1L == length(info) && is.na(info))
        info <- as.character(info)
    if (1L == length(geno) && is.na(geno))
        geno <- as.character(geno)
    if (asGRanges) {
       if (length(geno) == 0 && length(info) == 0)
           stop("when 'asGRanges=TRUE' one of 'geno' or 'info' ",
                "must be specified")
       if (length(geno) > 0 && length(info) > 0)
           stop("when 'asGRanges=TRUE' only one of 'geno' or 'info' ",
                "can be specified")
       if (length(geno) > 1)
           stop("when 'asGRanges=TRUE' only 1 element of 'geno' can ",
                "be specified")
    }
    new(class, which=which, info=info, geno=geno,
        trimEmpty=trimEmpty, asGRanges=asGRanges)
}

## ScanBcfParam

setMethod(ScanBcfParam, c(which="missing"),
    function(info=character(), geno=character(), trimEmpty=TRUE, which,
             asGRanges=FALSE, ...)
{
    which <- IRangesList()
    names(which) <- character()
    .ScanBcfParam(info, geno, trimEmpty, which, asGRanges, ...)
})

setMethod(ScanBcfParam, c(which="RangesList"), 
    function(info=character(), geno=character(), trimEmpty=TRUE, which,
             asGRanges=FALSE, ...)
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
             asGRanges=FALSE, ...)
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

