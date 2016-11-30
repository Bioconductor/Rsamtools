test_RsamtoolsFile_constructor <- function() {
    fl <- system.file("extdata", "example.gtf.gz", package="Rsamtools")
    checkException(TabixFile(c(fl, fl)), silent=TRUE)

    checkTrue(validObject(tbx <- TabixFile(fl, character())))
    checkIdentical(character(0), index(tbx, asNA=FALSE))

    checkTrue(validObject(tbx <- TabixFile(fl, NA)))
    checkIdentical(character(0), index(tbx, asNA=FALSE))
}

test_RsamtoolsFileList_constructor <- function() {
    fl <- system.file("extdata", "example.gtf.gz", package="Rsamtools")
    fls <- c(fl, fl)

    checkTrue(validObject(tbx <- TabixFileList(fls)))
    checkIdentical(setNames(fls, basename(fls)), path(tbx))
    checkIdentical(setNames(paste(fls, "tbi", sep="."), basename(fls)),
                   index(tbx))

    checkTrue(validObject(tbx <- TabixFileList(fls, character())))
    checkIdentical(setNames(fls, basename(fls)), path(tbx))
    checkIdentical(setNames(rep(NA_character_, 2), basename(fls)),
                   index(tbx))
    
    checkTrue(validObject(tbx <- TabixFileList(fls, NA)))
    checkIdentical(setNames(fls, basename(fls)), path(tbx))
    checkIdentical(setNames(rep(NA_character_, 2), basename(fls)),
                   index(tbx))
}
