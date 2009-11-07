Cigar <-
    function(cigars=factor(), ...)
{
    new("Cigar", .cigar=factor(cigars),  ...)
}

setMethod(initialize, "Cigar",
          function(.Object, ..., .cigar=cigars(.Object))
{
    callNextMethod(.Object, ..., .cigar=.cigar)
})

cigars <- function(x) slot(x, ".cigar")

.ucigars <- function(x) levels(cigars(x))

.cigartable <- function(x)
{
    .Call(.cigar_table, levels(cigars(x)))
}

setMethod(.validity, "Cigar", function(object) {
    msg <- NULL
    ucigar <- unique(cigars(object))
    if (!setequal(ucigar[!is.na(ucigar)], .ucigars(object)))
        msg <- c(msg,
                 "'.cigar' levels differ from elements")
    if (is.null(msg)) TRUE else msg
})

setMethod(length, "Cigar", function(x) length(cigars(x)))

setMethod(names, "Cigar", function(x) names(cigars(x)))

setMethod("[", c(x="Cigar", i="ANY", j="missing"),
          function(x, i, j, ..., drop=TRUE) 
{
    initialize(x, .cigar=factor(cigars(x)[i]))
})

setMethod("[[", c(x="Cigar", i="ANY", j="missing"),
          function(x, i, j, ...) as.character(cigars(x)[[i]]))

setMethod(show, "Cigar", function(object) {
    .show_classname(object)
    cat(sprintf("cigars (%d):", length(object)),
        as.character(IRanges:::selectSome(cigars(object))), "\n")
})
