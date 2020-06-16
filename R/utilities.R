.ppath <- function(tag, filepath)
{
    wd <- options('width')[[1]] - nchar(tag) - 6
    if(is.na(filepath))
        return(sprintf("%s: %s\n", tag, NA_character_))
    if (0L == length(filepath) || nchar(filepath) < wd)
        return(sprintf("%s: %s\n", tag, filepath))
    bname <- basename(filepath)
    wd1 <- wd - nchar(bname)
    dname <- substr(dirname(filepath), 1, wd1)
    sprintf("%s: %s...%s%s\n",
            tag, dname, .Platform$file.sep, bname)
}

.io_check_exists <-
    function(files)
{
    if (!length(files))
        stop("'files' is length(0)")
    idx <- !grepl("^(gs|aws|ftp|http|https)://", files) & !is.na(files)
    test <- file.exists(files[idx])
    if (!all(test)) {
        msg <- paste0(sQuote(files[idx][!test]), collapse = "\n  ")
        stop("file(s) do not exist:\n  ", msg)
    }
}

.show_classname <-
    function(x) cat("class: ", class(x), "\n", sep="")

.normalizePath <-
    function(path)
{
    if (is(path, "RsamtoolsFile")) {
        path <- path(path)
    } else {
        path <- as.character(path)
    }
    idx <- !grepl("^(ftp)|(http)://", path) & !is.na(path)
    ## expand ~/, but don't chase links (i.e., don't normalizePath())
    path[idx] <- path.expand(path[idx])
    path
}

.file.rename <-
    function(from, to)
{
    warn <- err <- NULL
    ok <- withCallingHandlers(tryCatch({
        file.rename(from, to) ||
            (file.copy(from, to) && file.remove(from))
    }, error=function(e) {
        err <<- append(err, conditionMessage(e))
        NULL
    }), warning=function(w) {
        warn <<- append(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    if (!ok) {
        msg <- "file.rename or file.copy/file.remove failed:\n  from: %s\n  to: %s\n  message(s): %s"
        stop(sprintf(msg, from, to, paste(c(warn, err), collapse="\n      ")))
    }
    ok
}

.uunlist <-
    function(x) unlist(x, use.names=FALSE)

setMethod(.asRegions, "IntegerRangesList", function(x) {
    list(as.character(space(x)), .uunlist(start(x)), .uunlist(end(x)))
})

setMethod(.asRegions, "GRanges", function(x) {
    list(as.character(seqnames(x)), start(x), end(x))
})

### All arguments must be parallel vectors (of length N).
### The arguments prefixed with 'x_' describe a vector 'x' of N alignments.
### The arguments prefixed with 'y_' describe a vector 'y' of N alignments.
### Performs "parallel pairing" of the N alignments in 'x' with the N
### alignments in 'y'.
.isValidHit <- function(x_flag, x_seqnames, x_start, x_mrnm, x_mpos,
                        y_flag, y_seqnames, y_start, y_mrnm, y_mpos)
{
    .Call(.p_pairing, NULL, x_flag, x_seqnames, x_start, x_mrnm, x_mpos,
                      NULL, y_flag, y_seqnames, y_start, y_mrnm, y_mpos)
}

### 'x_flag', 'x_seqnames', 'x_start', 'x_mrnm', 'x_mpos': parallel vectors
###     (of length N) describing N alignments. The alignments are assumed to
###     be already grouped by QNAME.
### 'group.sizes': vector of non-negative integers which sum to N.
###     If 'x_qname' was a vector of length N parallel to the 'x_*' arguments
###      and containing the QNAME field, then 'group.sizes' would be
###     'runLength(Rle(x_qname))'.
### Returns an integer vector of length N parallel to the 'x_*' arguments.
### Alignments with more than 1 possible mate are assigned a zero.
### Those with exactly 1 mate that has itself more than 1 mate are assigned
### a negative value (the opposite of the index of the mate).
.findMateWithinGroups <- function(group.sizes,
                                  x_flag, x_seqnames,
                                  x_start, x_mrnm, x_mpos)
{
    .Call(.find_mate_within_groups, group.sizes,
                                    x_flag, x_seqnames,
                                    x_start, x_mrnm, x_mpos)
}

