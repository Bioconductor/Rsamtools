\name{BcfInput}
\Rdversion{1.1}

\alias{scanBcfHeader}
\alias{scanBcfHeader,character-method}
\alias{scanBcf}
\alias{scanBcf,character-method}
\alias{asBcf}
\alias{asBcf,character-method}
\alias{indexBcf}
\alias{indexBcf,character-method}


\title{

  Operations on `BCF' files.

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

}

\arguments{

  \item{file}{For \code{scanBcf} and \code{scanBcfHeader}, the
    character() file name of the \sQuote{BCF} file to
    be processed, or an instance of class \code{\link{BcfFile}}.} 

  \item{index}{The character() file name(s) of the `BCF' index to be
    processed.}

  \item{dictionary}{a character vector of the unique \dQuote{CHROM}
    names in the VCF file.}

  \item{destination}{The character(1) file name of the location where
    the BCF output file will be created. For \code{asBcf} this is
    without the \dQuote{.bcf} file suffix.}

  \item{param}{A instance of \code{\linkS4class{ScanBcfParam}}
    influencing which records are parsed and the \sQuote{INFO} and 
    \sQuote{GENO} information returned.}

  \item{...}{Additional arguments, e.g., for
    \code{scanBcfHeader,character-method}, \code{mode} of
    \code{\link{BcfFile}}.}

  \item{overwrite}{A logical(1) indicating whether the destination
    can be over-written if it already exists.}

  \item{indexDestination}{A logical(1) indicating whether the created
    destination file should also be indexed.}

}

\details{

  \code{bcf*} functions are restricted to the GENO fields supported by 
  \sQuote{bcftools} (see documentation at the url below). The argument 
  \code{param} allows portions of the file to be input, but requires 
  that the file be BCF or bgzip'd and indexed as a 
  \code{\linkS4class{TabixFile}}. For similar functions operating on VCF
  files see ?\code{scanVcf} in the \code{VariantAnnotation} package.

}

\value{

  \code{scanBcfHeader} returns a list, with one element for each file 
  named in \code{file}. Each element of the list is itself a list containing 
  three elements. The \code{reference} element is a character() vector with 
  names of reference sequences. The \code{sample} element is a character() 
  vector of names of samples. The \code{header} element is a character() 
  vector of the header lines (preceeded by \dQuote{##}) present in the VCF 
  file.

  \code{scanBcf} returns a list, with one element per file. Each list has 9 
  elements, corresponding to the columns of the VCF specification: \code{CHROM}, 
  \code{POS}, \code{ID}, \code{REF}, \code{ALT}\code{QUAL}, \code{FILTER}, 
  \code{INFO}, \code{FORMAT}, \code{GENO}.

  The \code{GENO} element is itself a list, with elements corresponding
  to fields supported by \sQuote{bcftools} (see documentation at the url below).

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
