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
    checkIdentical(TRUE, file.remove(d))
}
