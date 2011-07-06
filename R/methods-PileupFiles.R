PileupFiles <-
    function(files, ..., param=PileupParam())
{
    bfl <- BamFileList(files, ...)
    new("PileupFiles", bfl, param=param)
}

plpFiles <- function(object) as(object, "BamFileList")

plpParam <- function(object) object@param

setMethod(applyPileups, c("PileupFiles", "PileupParam"),
    function(files, FUN, ..., param)
{
    FUN <- match.fun(FUN)
    ok <- isOpen(files)
    if (!all(ok))
        if (any(ok))
            stop("all(isOpen(<PileupFiles>))' is not 'TRUE'")
        else {
            open(files)
            on.exit(close(files))
        }
    tryCatch({
        param <- as(param, "list")
        extptr <- lapply(files, .extptr)
        space <-
            if (0L != length(param[["which"]])) .asSpace(param[["which"]])
            else NULL
        what <- logical(2);
        param[["what"]] <- c("seq", "qual") %in% param[["what"]]
        .Call(.apply_pileups, extptr, names(files), space, param, FUN)
    }, error=function(err) {
        stop("applyPileups: ", conditionMessage(err), call.=FALSE)
    })
})

setMethod(applyPileups, c("PileupFiles", "missing"),
    function(files, FUN, ..., param)
{
    applyPileups(files, FUN, ..., param=plpParam(files))
})

setMethod(show, "PileupFiles", function(object) {
    cat("class:", class(object), "\n")
    nms <- names(object)
    txt <- paste(IRanges:::selectSome(nms, 3L), collapse=", ")
    cat(sprintf("names: %s (%d total)\n", txt, length(nms)))
    fls <- sapply(object, function(x) basename(path(x)))
    txt <- paste(IRanges:::selectSome(fls, 3L), collapse=", ")
    cat(sprintf("plpFiles: %s (%d total)\n", txt, length(fls)))
    cat("plpParam: class", class(plpParam(object)), "\n")
})
