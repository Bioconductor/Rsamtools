setMethod(.validity, "PileupFiles", function(object) {
    msg <- NULL
    ok <- sapply(object$.files, is, "RsamtoolsFile")
    if (!all(ok))
        msg <- c(msg, "all elements must extend 'RsamtoolsFile'")
    if (is.null(msg)) TRUE else msg
})

setMethod(PileupFiles, "character",
    function(..., param=PileupParam())
{
    fls <- lapply(..1, function(x) open(BamFile(x)))
    res <- .PileupFiles$new(files=fls, param=param)
    validObject(res)
    res
})

setMethod(PileupFiles, "BamFile",
    function(..., param=PileupParam())
{
    fls <- lapply(list(...), function(x) {
        if (!isOpen(x)) open(x) else x
    })
    .PileupFiles$new(files=fls, param=param)
})

open.PileupFiles <-
    function(con, ...)
{
    for (fl in con$files)
        open(fl)
    invisible(con)
}

close.PileupFiles <-
    function(con, ...)
{
    for (fl in con$files)
        close(fl)
    invisible(con)
}

setMethod(isOpen, "PileupFiles",
    function(con, rw="")
{
    if (!missing(rw) && rw != "read")
        stop("'rw' must be 'read'")
    sapply(con$files, isOpen)
})

plpFiles <- function(object) object$files

plpParam <- function(object) object$param

setMethod(applyPileups, c("PileupFiles", "PileupParam"),
    function(files, FUN, ..., param)
{
    FUN <- match.fun(FUN)
    tryCatch({
        param <- as(param, "list")
        files <- lapply(files$files, .extptr)
        space <-
            if (0L != length(param[["which"]])) .asSpace(param[["which"]])
            else NULL
        what <- logical(2);
        param[["what"]] <- c("seq", "qual") %in% param[["what"]]
        .Call(.apply_pileups, files, space, param, FUN)
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
    fls <- sapply(plpFiles(object), function(x) basename(path(x)))
    txt <- paste(selectSome(fls, 3L), collapse=", ")
    cat(sprintf("plpFiles: %s (%d total)\n", txt, length(fls)))
    cat("plpParam: class", class(object$param), "\n")
})
