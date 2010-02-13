test_Alignments0_constructor <- function()
{
    checkTrue(validObject(Alignments0()))
    checkTrue(validObject(Alignments0(rname=factor("A"),
                                      strand=strand("-"),
                                      pos=1L, cigar="1M")))
    checkTrue(validObject(Alignments0(rname=factor("A"),
                                      strand=strand(NA_character_),
                                      pos=1L, cigar="1M")))
}

test_Alignments1_constructor <- function()
{
    checkTrue(validObject(Alignments1()))
    checkTrue(validObject(Alignments1(rname=factor("A"),
                                      strand=strand("-"),
                                      pos=1L, cigar="1M")))
    checkTrue(validObject(Alignments1(rname=factor("A"),
                                      strand=strand(NA_character_),
                                      pos=1L, cigar="1M")))
}
