setMethod(.validity, "PileupParam", function(object) {
    msg <- NULL
    flag <- object@flag
    if (2L != length(flag) ||
        any(c("keep0", "keep1") != names(flag)))
        msg <- c(msg, "'flag' not from scanBamFlag()")
    len1elts <- c("minBaseQuality", "minMapQuality", "minDepth",
                  "maxDepth", "yieldSize", "yieldBy", "yieldAll")
    ok <- 1L == sapply(len1elts,
            function(x, obj) length(slot(obj, x)),
            object)
    if (!all(ok))
        msg <- c(msg,
                 sprintf("'%s' must be length 1",
                         paste(len1elts[!ok], collapse="' '")))
    what <- eval(formals(PileupParam)[["what"]])
    ok <- object@what %in% what
    if (!all(ok))
        msg <- c(msg, sprintf("'what' must be in '%s'",
                              paste(what, collapse="' '")))
    if (is.null(msg)) TRUE else msg
})

PileupParam <-
    function(flag=scanBamFlag(),
             minBaseQuality=13L,
             minMapQuality=0L,
             minDepth=0L,
             maxDepth=250L,
             yieldSize=1L,
             yieldBy=c("range", "position"),
             yieldAll=FALSE,
             which=GRanges(),
             what=c("seq", "qual"))
{
    yieldBy <- match.arg(yieldBy)
    if ("range" == yieldBy && yieldSize != 1)
        stop("'yieldSize' must equal 1 when 'yieldBy=\"range\"'")
    new("PileupParam", flag=flag,
        minBaseQuality=as.integer(minBaseQuality),
        minMapQuality=as.integer(minMapQuality),
        minDepth=as.integer(minDepth),
        maxDepth=as.integer(maxDepth),
        yieldSize=as.integer(yieldSize),
        yieldBy=yieldBy,
        yieldAll=as.logical(yieldAll),
        which=which, what=what)
}

setAs("PileupParam", "list", function(from) {
    nms <- slotNames(class(from))
    res <- lapply(nms, slot, object=from)
    names(res) <- nms
    res
})

plpFlag <- function(object) slot(object, "flag")
"plpFlag<-" <- function(object, value)
{
    slot(object, "flag") <- value
    validObject(object)
    object
}

plpMinBaseQuality <-function(object) slot(object, "minBaseQuality")
"plpMinBaseQuality<-" <- function(object, value)
{
    slot(object, "minBaseQuality") <- as.integer(value)
    validObject(object)
    object
}

plpMinMapQuality <-
    function(object) slot(object, "minMapQuality")
"plpMinMapQuality<-" <- function(object, value)
{
    slot(object, "minMapQuality") <- as.integer(value)
    validObject(object)
    object
}

plpMinDepth <- function(object) slot(object, "minDepth")
"plpMinDepth<-" <- function(object, value)
{
    slot(object, "minDepth") <- as.integer(value)
    validObject(object)
    object
}

plpMaxDepth <- function(object) slot(object, "maxDepth")
"plpMaxDepth<-" <- function(object, value)
{
    slot(object, "maxDepth") <- as.integer(value)
    validObject(object)
    object
}

plpYieldSize <- function(object) slot(object, "yieldSize")
"plpYieldSize<-" <- function(object, value)
{
    slot(object, "yieldSize") <- as.integer(value)
    validObject(object)
    object
}

plpYieldBy <- function(object) slot(object, "yieldBy")
"plpYieldBy<-" <- function(object, value)
{
    slot(object, "yieldBy") <- as.character(value)
    validObject(object)
    object
}

plpYieldAll <- function(object) slot(object, "yieldAll")
"plpYieldAll<-" <- function(object, value)
{
    slot(object, "yieldAll") <- as.logical(value)
    validObject(object)
    object
}

plpWhich <- function(object) slot(object, "which")
"plpWhich<-" <- function(object, value)
{
    slot(object, "which") <- as(value, "GRanges")
    validObject(object)
    object
}

plpWhat <- function(object) slot(object, "what")
"plpWhat<-" <- function(object, value)
{
    slot(object, "what") <- as.character(value)
    validObject(object)
    object
}
    
setMethod(show, "PileupParam", function(object) {
    cat("class:", class(object), "\n")
    cat("plpFlag:",
        sprintf("%s=%s", names(object@flag), object@flag),
        "\n")
    len1elts <- c("minBaseQuality", "minMapQuality", "minDepth",
                  "maxDepth", "yieldSize", "yieldBy", "yieldAll")
    for (elt in len1elts)
        cat(sprintf("%s: %s",
                    sub("([[:alpha:]])", "plp\\U\\1", elt, perl=TRUE),
                    slot(object, elt)), "\n")
    cat("plpWhat: '", paste(object@what, collapse="' '"), "'\n", sep="")
    cat(sprintf("plpWhich: %s (length %d)\n", class(object@which),
                length(object@which)))
})
