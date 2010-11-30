\name{BcfInput}
\Rdversion{1.1}

\alias{scanBcfHeader}
\alias{scanBcfHeader,character-method}
\alias{scanBcf}
\alias{scanBcf,character-method}

\title{

  Operations on `VCF' or `BCF' (variant call) files.

}
\description{
  Scan headers of variant call files in text or binary format.
}
\usage{

scanBcfHeader(file, ...)
\S4method{scanBcfHeader}{character}(file, ...)

scanBcf(file, ...)
\S4method{scanBcf}{character}(file, index = file, ..., param=ScanBcfParam())

}

\arguments{

  \item{file}{The character() file name(s) of the `VCF' or `BCF' file to
    be processed, or an instance of class \code{\link{BcfFile}}.}

  \item{index}{The character() file name(s) of the `BCF' index to be
    processed.}

  \item{param}{A instance of \code{\linkS4class{ScanBcfParam}}
    influencing which records are parsed and the \sQuote{INFO} and
    \sQuote{GENO} information returned.}

  \item{...}{Additional arguments, e.g., for
    \code{scanBcfHeader,character-method}, \code{mode} of
    \code{\link{BcfFile}}.}

}

\value{

  \code{scanBcfHeader} returns a list, with one element for each file
  named in \code{files}. Each element of the list is itself a list
  containing three element. The \code{reference} element is a
  character() vector with names of reference sequences.  The
  \code{sample} element is a character() vector of names of samples. The
  \code{header} element is a character() vector of the header lines
  (preceeded by \dQuote{##}) present in the VCF file.

  \code{scanBcf} returns a list, with one element per file. Each list
  has 9 elements, corresponding to the columns of the VCF specification:
  \code{CHROM}, \code{POS}, \code{ID}, \code{REF},
  \code{ALT}\code{QUAL}, \code{FILTER}, \code{INFO}, \code{FORMAT},
  \code{GENO}. The \code{GENO} element is itself a list, with elements
  corresponding to those supported by \sQuote{bcftools} (see
  documentation at the url below).

}

\references{


  \url{http://vcftools.sourceforge.net/specs.html} outlines the VCF
  specification.

  \url{http://samtools.sourceforge.net/mpileup.shtml} contains
  information on the portion of the specification implemented by
  \code{bcftools}.

  \url{http://samtools.sourceforge.net/} provides information on
  \code{samtools}.

}

\author{

  Martin Morgan <mtmorgan@fhcrc.org>.

}

\examples{
fl <- system.file("extdata", "ex1.bcf", package="Rsamtools")
scanBcfHeader(fl)
bcf <- scanBcf(fl)
## value: list-of-lists
str(bcf[1:8])
names(bcf[["GENO"]])
str(head(bcf[["GENO"]][["PL"]]))
example(BcfFile)
}

\keyword{ manip }
