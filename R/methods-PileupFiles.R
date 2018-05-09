PileupFiles <-
    function(files, ..., param=ApplyPileupsParam())
{
    bfl <- BamFileList(files, ...)
    new("PileupFiles", bfl, param=param)
}

plpFiles <- function(object) as(object, "BamFileList")

plpParam <- function(object) object@param

setMethod(applyPileups, c("PileupFiles", "ApplyPileupsParam"),
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
    lvls <- lapply(files, seqlevels)
    for (i in seq_along(files)[-1])
        if (!identical(lvls[[i]], lvls[[1]])) {
            msg <- sprintf("applyPileups 'seqlevels' must be identical();
                failed when comparing %s with %s",
                sQuote(basename(path(files)[1])),
                sQuote(basename(path(files)[i])))
            stop(paste(strwrap(msg, exdent=4), collapse="\n"))
        }
    tryCatch({
        param <- as(param, "list")
        extptr <- lapply(files, .extptr)
        regions <-
            if (0L != length(param[["which"]])) .asRegions(param[["which"]])
            else NULL
        param[["what"]] <- c("seq", "qual") %in% param[["what"]]
        .Call(.apply_pileups, extptr, names(files), regions, param, FUN)
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
    txt <- paste(S4Vectors:::selectSome(nms, 3L), collapse=", ")
    cat(sprintf("names: %s (%d total)\n", txt, length(nms)))
    fls <- sapply(object, function(x) basename(path(x)))
    txt <- paste(S4Vectors:::selectSome(fls, 3L), collapse=", ")
    cat(sprintf("plpFiles: %s (%d total)\n", txt, length(fls)))
    cat("plpParam: class", class(plpParam(object)), "\n")
})
