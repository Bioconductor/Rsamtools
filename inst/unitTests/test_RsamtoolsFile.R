test_RsamtoolsFile_constructor <- function() {
    fl <- system.file("extdata", "example.gtf.gz", package="Rsamtools")
    checkException(TabixFile(c(fl, fl)), silent=TRUE)

    checkTrue(validObject(tbx <- TabixFile(fl, character())))
    checkIdentical(index(tbx), character())

    checkTrue(validObject(tbx <- TabixFile(fl, NA)))
    checkIdentical(index(tbx), character())
}

test_RsamtoolsFileList_constructor <- function() {
    fl <- system.file("extdata", "example.gtf.gz", package="Rsamtools")
    fls <- c(fl, fl)

    checkTrue(validObject(tbx <- TabixFileList(fls)))
    checkIdentical(setNames(fls, basename(fls)), path(tbx))
    checkIdentical(setNames(paste(fls, "tbi", sep="."), basename(fls)),
                   vapply(as.list(tbx), index, character(1)))

    checkTrue(validObject(tbx <- TabixFileList(fls, character())))
    checkIdentical(setNames(fls, basename(fls)), path(tbx))
    checkIdentical(setNames(rep(list(character(0)), 2), basename(fls)),
                   sapply(as.list(tbx), index))
    
    checkTrue(validObject(tbx <- TabixFileList(fls, NA)))
    checkIdentical(setNames(fls, basename(fls)), path(tbx))
    checkIdentical(setNames(rep(list(character(0)), 2), basename(fls)),
                   sapply(as.list(tbx), index))
}
