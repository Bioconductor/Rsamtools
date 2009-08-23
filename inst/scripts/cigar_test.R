library(Rsamtools)
library(RUnit)

res <- local({
    fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
    scanBam(fl)[[1]]
})

test_cigar_run_count <- function() {
    ucig <- unique(as.character(cigars(res[["cigar"]])))
    lens <- .Call(Rsamtools:::.cigar_run_count, ucig)
    checkEquals(length(ucig), length(lens))
    checkEquals(56L, sum(lens))
    checkEquals(1L, sum(lens == 0))
}

test_cigar_table <- function() {
    ucig <- unique(as.character(cigars(res[["cigar"]])))
    df <- .Call(Rsamtools:::.cigar_table, ucig)
    checkEquals(c("cigar", "element", "length", "value"), names(df))
    checkEquals(ucig, unique(df[["cigar"]]))

    lens <- .Call(Rsamtools:::.cigar_run_count, ucig)
    checkEquals(sum(lens) + sum(lens==0), nrow(df))
}
