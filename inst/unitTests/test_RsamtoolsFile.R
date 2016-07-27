test_RsamtoolsFile_constructor <- function() {
    fl <- system.file("extdata", "example.gtf.gz", package="Rsamtools")
    checkException(TabixFile(c(fl, fl)), silent=TRUE)
}
