test_sortBam <- function() {
    fl0 <- system.file("extdata", "ex1.bam", package="Rsamtools")
    src <- system.file("unitTests", "cases", package="Rsamtools")
    fl <- file.path(src, "ex1_unsort.bam")
    ofl <- tempfile()
    checkTrue(file.create(ofl))
    on.exit(unlink(ofl))
    sorted <- sortBam(fl, ofl)
    exp <- scanBam(fl0)[[1]]
    obs <- scanBam(sorted)[[1]]
    checkIdentical(exp[["rname"]], obs[["rname"]])
    checkIdentical(Filter(Negate(is.na), exp[["pos"]]),
                   Filter(Negate(is.na), obs[["pos"]]))
}

test_sortBam_not_BAM_input <- function() {
    fl0 <- system.file("extdata", "ex1.sam", package="Rsamtools")
    checkException(sortBam(fl0, tempfile()), silent=TRUE)
}

test_sortBam_byTag <- function() {
  src <- system.file("unitTests", "cases", package="Rsamtools")
  fl <- file.path(src, "ex1_unsort.bam")
  ofl <- tempfile()
  checkTrue(file.create(ofl))
  on.exit(unlink(ofl))

  # sort by integer Tag Aq
  sorted <- sortBam(fl, ofl, byTag = "Aq")
  obs <- scanBam(sorted, param = ScanBamParam(tag = "Aq"))[[1]]
  tVal <- obs$tag$Aq

  # reads without Aq tag are first records in sorted bam
  checkIdentical(which(is.na(tVal)), 1L:36L)
  validTags <- tVal[!is.na(tVal)]
  checkIdentical(validTags, sort(validTags))

  checkException(sortBam(fl, ofl, byTag = 1), silent=TRUE)
  checkException(sortBam(fl, ofl, byTag = c("bogus", "input")), silent=TRUE)
}

test_sortBam_nThreads <- function() {
  src <- system.file("unitTests", "cases", package="Rsamtools")
  fl <- file.path(src, "ex1_unsort.bam")
  checkException(sortBam(fl, tempfile(), nThreads = 0), silent=TRUE)
  checkException(sortBam(fl, tempfile(), nThreads = c(0, 1)), silent=TRUE)
}
