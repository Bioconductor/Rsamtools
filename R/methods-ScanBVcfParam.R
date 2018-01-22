.ScanBcfParam <-
    function(fixed=character(), info=character(), geno=character(), 
             samples=character(), trimEmpty=TRUE, which, 
             class="ScanBcfParam")
{
    if (1L == length(fixed) && is.na(fixed))
        fixed <- as.character(fixed)
    if (1L == length(info) && is.na(info))
        info <- as.character(info)
    if (1L == length(geno) && is.na(geno))
        geno <- as.character(geno)
    if (1L == length(samples) && is.na(samples))
        samples <- as.character(samples)
    new(class, which=which, fixed=fixed, info=info, geno=geno,
        samples=samples, trimEmpty=trimEmpty)
}

## ScanBcfParam

setMethod(ScanBcfParam, c(which="missing"),
    function(fixed=character(), info=character(), geno=character(), 
             samples=character(), trimEmpty=TRUE, which, ...)
{
    which <- IRangesList()
    names(which) <- character()
    .ScanBcfParam(fixed, info, geno, samples, trimEmpty, which, ...)
})

setMethod(ScanBcfParam, c(which="GRangesList"), 
    function(fixed=character(), info=character(), geno=character(), 
             samples=character(), trimEmpty=TRUE, which, ...)
{
    .ScanBcfParam(fixed, info, geno, samples, trimEmpty, 
                  which=ranges(which), ...)
})

setMethod(ScanBcfParam, c(which="IntegerRangesList"), 
    function(fixed=character(), info=character(), geno=character(), 
             samples=character(), trimEmpty=TRUE, which, ...)
{
    .ScanBcfParam(fixed, info, geno, samples, trimEmpty, which, ...)
})

setMethod(ScanBcfParam, c(which="GRanges"),
    function(fixed=character(), info=character(), geno=character(), 
             samples=character(), trimEmpty=TRUE, which, ...)
{
    which <- split(ranges(which), seqnames(which))
    .ScanBcfParam(fixed, info, geno, samples, trimEmpty, which, ...)
})

## accessors

bcfFixed <- function(object) slot(object, "fixed")
bcfInfo <- function(object) slot(object, "info")
bcfGeno <- function(object) slot(object, "geno")
bcfSamples <- function(object) slot(object, "samples")
bcfTrimEmpty <- function(object) slot(object, "trimEmpty")
bcfWhich <- function(object) slot(object, "which")

.some <- S4Vectors:::selectSome
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
    cat(.clslbl("Info:"), .some(bcfInfo(object)), "\n")
    cat(.clslbl("Geno:"), .some(bcfGeno(object)), "\n")
    cat(.clslbl("Samples:"), .some(bcfSamples(object)), "\n")
})

