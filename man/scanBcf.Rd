\name{VcfInput}
\Rdversion{1.1}

\alias{scanBcfHeader}
\alias{scanBcfHeader,character-method}
\alias{scanBcf}
\alias{scanBcf,character-method}
\alias{asBcf}
\alias{asBcf,character-method}
\alias{indexBcf}
\alias{indexBcf,character-method}

\alias{scanVcfHeader}
\alias{scanVcfHeader,character-method}
\alias{scanVcf}
\alias{scanVcf,character,ANY-method}
\alias{scanVcf,character,missing-method}
\alias{scanVcf,connection,missing-method}

\title{

  Operations on `VCF' or `BCF' (variant call) files.

}
\description{
  Import, coerce, or index variant call files in text or binary format.
}
\usage{

scanBcfHeader(file, ...)
\S4method{scanBcfHeader}{character}(file, ...)

scanBcf(file, ...)
\S4method{scanBcf}{character}(file, index = file, ..., param=ScanBcfParam())

asBcf(file, dictionary, destination, ...,
      overwrite=FALSE, indexDestination=TRUE)
\S4method{asBcf}{character}(file, dictionary, destination, ...,
      overwrite=FALSE, indexDestination=TRUE)

indexBcf(file, ...)
\S4method{indexBcf}{character}(file, ...)

scanVcfHeader(file, ...)
\S4method{scanVcfHeader}{character}(file, ...)

scanVcf(file, ..., param)
\S4method{scanVcf}{character,ANY}(file, ..., param)
\S4method{scanVcf}{character,missing}(file, ..., param)
\S4method{scanVcf}{connection,missing}(file, ..., param)
}

\arguments{

  \item{file}{For \code{scanBcf} and \code{scanBcfHeader}, the
    character() file name of the \sQuote{VCF} or \sQuote{BCF} file to
    be processed, or an instance of class \code{\link{BcfFile}}. For
    \code{scanVcf} and \code{scanVcfHeader}, the character() file name,
    \code{\link{TabixFile}}, or class \code{connection} ( \code{file()}
    or \code{bgzip()}) of the \sQuote{VCF} file to be processed.}

  \item{index}{The character() file name(s) of the `BCF' index to be
    processed.}

  \item{dictionary}{a character vector of the unique \dQuote{CHROM}
    names in the VCF file.}

  \item{destination}{The character(1) file name of the location where
    the BCF output file will be created. For \code{asBcf} this is
    without the \dQuote{.bcf} file suffix.}

  \item{param}{A instance of \code{\linkS4class{ScanBcfParam}} or
    \code{\linkS4class{ScanVcfParam}} influencing which records are
    parsed and the \sQuote{INFO} and \sQuote{GENO} information
    returned.}

  \item{...}{Additional arguments, e.g., for
    \code{scanBcfHeader,character-method}, \code{mode} of
    \code{\link{BcfFile}}.}

  \item{overwrite}{A logical(1) indicating whether the destination
    can be over-written if it already exists.}

  \item{indexDestination}{A logical(1) indicating whether the created
    destination file should also be indexed.}

}

\details{

  Most users will use the \code{vcf*} functions; \code{bcf*} are
  restricted to the GENO fields supported by \sQuote{bcftools} (see
  documentation at the url below). The argument \code{param} allows
  protions of the file to be input, but requires that the file be BCF or
  bgzip'd and indexed as a \code{\linkS4class{TabixFile}}.

  \code{scanVcf} with \code{param="missing"} and \code{file="character"}
  or \code{file="connection"} parse the entire file. With
  \code{file="connection"}, an argument \code{n} indicates the number of
  lines of the VCF file to input; a connection open at the begining of
  the call is open and incremented by \code{n} lines at the end of the
  call, providing a convenient way to stream through large VCF files.

}

\value{

  \code{scanVcfHeader} / \code{scanBcfHeader} returns a list, with one
  element for each file named in \code{file}. Each element of the list
  is itself a list containing three element. The \code{reference}
  element is a character() vector with names of reference sequences.
  The \code{sample} element is a character() vector of names of
  samples. The \code{header} element is a character() vector of the
  header lines (preceeded by \dQuote{##}) present in the VCF file.

  \code{scanVcf} / \code{scanBcf} returns a list, with one element per
  file. Each list has 9 elements, corresponding to the columns of the
  VCF specification: \code{CHROM}, \code{POS}, \code{ID}, \code{REF},
  \code{ALT}\code{QUAL}, \code{FILTER}, \code{INFO}, \code{FORMAT},
  \code{GENO}.

  The \code{GENO} element is itself a list, with elements corresponding
  to those defined in the VCF file header. For \code{scanVcf}, elements
  of GENO are returned as a matrix of records x samples; if the
  description of the element in the file header indicated multiplicity
  other than 1 (e.g., variable number for \dQuote{A}, \dQuote{G}, or
  \dQuote{.}), then each entry in the matrix is a character string with
  sub-entries comma-delimited.

  \code{asBcf} creates a binary BCF file from a text VCF file.

  \code{indexBcf} creates an index into the BCF file.

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

\seealso{

  \code{\link{BcfFile}}, \code{\link{TabixFile}}

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
