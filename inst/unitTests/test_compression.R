test_bgzip_openclose <- function()
{
    ## trying to determine that file handle has been cleaned up
    checkIdentical(TRUE, dir.create(d <- tempfile()))
    fin <- file.path(d, "in")
    fout <- file.path(d, "out")
    writeLines("123", con=fin)
    bgzip(fin, fout)
    checkIdentical(TRUE, file.remove(fin))
    checkIdentical(TRUE, file.remove(fout))
    checkIdentical(0L, unlink(d, recursive=TRUE))
}

test_razip_small_files <- function()
{
    src <- system.file("extdata", "ce2dict1.fa", package="Rsamtools")
    file.copy(src, dest <- tempfile())
    checkIdentical(readLines(src), readLines(dest))
}

test_razip_gzcompressed <- function()
{
    file.copy(system.file("extdata", "ce2dict1.fa", package="Rsamtools"),
              src <- tempfile())
    con <- gzfile(gzdest <- tempfile(), "wb")
    writeLines(readLines(src), con, sep="\n")
    close(con)
    rzdest <- razip(gzdest)
    checkIdentical(src, indexFa(src))
    checkIdentical(gzdest, indexFa(gzdest))
    checkIdentical(rzdest, indexFa(rzdest))
    src <- scanFa(src, param=as(seqinfo(FaFile(rzdest)), "GRanges"))
    rz <- scanFa(rzdest, param=as(seqinfo(FaFile(rzdest)), "GRanges"))
    gz <- scanFa(rzdest, param=as(seqinfo(FaFile(rzdest)), "GRanges"))
    checkIdentical(as.character(src), as.character(rz))
    checkIdentical(as.character(src), as.character(gz))
}
