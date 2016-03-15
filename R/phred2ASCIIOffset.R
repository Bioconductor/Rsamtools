.ascii_offset <- function()
    setNames(33:126 - 33L, strsplit(rawToChar(as.raw(33:126)), "")[[1]])

.phred2ascii_int <-
    function(x, scheme)
{
    ## See https://en.wikipedia.org/wiki/FASTQ_format#Encoding
    ascii <- .ascii_offset()
    switch(scheme, "Illumina 1.8+" = {
        ## L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
        stopifnot(all(x >= 0), all(x <= 41))
        ascii[x + 1L]
    }, "Sanger" = {
        ## S - Sanger        Phred+33,  raw reads typically (0, 40)
        stopifnot(all(x >= 0), all(x <= 40))
        ascii[x + 1L]
    }, "Solexa" = {
        ## X - Solexa        Solexa+64, raw reads typically (-5, 40)
        stopifnot(all(x >= -5), all(x <= 40))
        ascii[x + 32L]
    }, "Illumina 1.3+" = {
        ## I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
        stopifnot(all(x >= 0), all(x <= 40))
        ascii[x + 32L]
    }, "Illumina 1.5+" = {
        ## J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
        ##     with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold) 
        ##     (Note: See discussion above).
        stopifnot(all(x >= 3), all(x <= 40))
        ascii[x + 32L]
    }, default = stop("unknown scheme '", scheme, "'"))
}

.phred2ascii_char <-
    function(x)
{
    ascii <- .ascii_offset()
    stopifnot(all(x %in% names(ascii)))
    ascii[x]
}

phred2ASCIIOffset <-
    function(phred=integer(),scheme= c("Illumina 1.8+", "Sanger", "Solexa",
                                       "Illumina 1.3+", "Illumina 1.5+"))
{
    if (is.numeric(phred)) {
        stopifnot(missing(scheme) || (length(scheme) == 1L), !anyNA(phred))
        scheme <- match.arg(scheme)
        phred <- as.integer(phred)
        .phred2ascii_int(phred, scheme)
    } else if (is.character(phred)) {
        if (!missing(scheme))
            message("'scheme' ignored; does not influence ASCII offset")
        if (length(phred) == 1L && nchar(phred) > 1L)
            phred <- strsplit(phred, "")[[1]]
        stopifnot(all(nchar(phred) == 1L))
        .phred2ascii_char(phred)
    } else
        stop("'phred' must be numeric (coerced to integer) or character")
}
