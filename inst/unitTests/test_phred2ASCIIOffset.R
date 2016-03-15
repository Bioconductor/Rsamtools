test_phred2ASCIIOffset <- function() {
    schemes <- eval(formals(phred2ASCIIOffset)$scheme)

    exp <- structure(integer(), names=character())
    checkIdentical(exp, phred2ASCIIOffset())
    checkIdentical(exp, phred2ASCIIOffset(character()))
    checkIdentical(exp, phred2ASCIIOffset(integer()))
    for (scheme in schemes) {
        checkIdentical(exp, phred2ASCIIOffset(scheme=scheme))
        checkIdentical(exp, phred2ASCIIOffset(integer(), scheme=scheme))
    }

    offsets <- Rsamtools:::.ascii_offset()
    checkIdentical(offsets[1L + 0:40], phred2ASCIIOffset(0:40, "Sanger"))
    checkIdentical(offsets[32L + (-5):40], phred2ASCIIOffset((-5):40, "Solexa"))
    checkIdentical(offsets[32L + 0:40],
                   phred2ASCIIOffset(0:40, "Illumina 1.3+"))
    checkIdentical(offsets[32L + 3:40],
                   phred2ASCIIOffset(3:40, "Illumina 1.5+"))
    checkIdentical(offsets[1L + 0:41], phred2ASCIIOffset(0:41, "Illumina 1.8+"))

    checkIdentical(offsets,
                   phred2ASCIIOffset(paste(names(offsets), collapse="")))

    checkException(phred2ASCIIOffset(50))
    checkException(phred2ASCIIOffset(-1))
    checkException(phred2ASCIIOffset(""))
    checkException(phred2ASCIIOffset(c("A", "AA")))
}
