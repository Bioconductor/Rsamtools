\name{indexTabix}
\Rdversion{1.1}

\alias{indexTabix}

\title{

  Compress and index tabix-compatible files.

}
\description{
  
  Index (with \code{indexTabix}) files that have been sorted into
  ascending sequence, start and end position ordering.

}
\usage{

indexTabix(file,
           format=c("gff", "bed", "sam", "vcf", "vcf4", "psltbl"),
           seq=integer(), start=integer(), end=integer(),
           skip=0L, comment="#", zeroBased=FALSE, ...)

}

\arguments{

  \item{file}{A characater(1) path to a sorted, bgzip-compressed file.}

  \item{format}{The format of the data in the compressed file. A
    characater(1) matching one of the types named in the function
    signature.}

  \item{seq}{If \code{format} is missing, then \code{seq} indicates the
    column in which the \sQuote{sequence} identifier (e.g., chrq) is to
    be found.}

  \item{start}{If \code{format} is missing, \code{start} indicates the
    column containing the start coordinate of the feature to be
    indexed.}

  \item{end}{If \code{format} is missing, \code{end} indicates the
    column containing the ending coordinate of the feature to be
    indexed.}

  \item{skip}{The number of lines to be skipped at the beginning of the
    file.}

  \item{comment}{A single character which, when present as the first
    character in a line, indicates that the line is to be omitted. from
    indexing.}

  \item{zeroBased}{A logical(1) indicating whether coordinats in the
    file are zero-based.}

  \item{...}{Additional arguments.}

}

\value{

  The return value of \code{indexTabix} is an updated instance of
  \code{file} reflecting the newly-created index file.

}

\references{

  \url{http://samtools.sourceforge.net/tabix.shtml}

}

\author{

  Martin Morgan <mtmorgan@fhcrc.org>.

}

\examples{
from <- system.file("extdata", "ex1.sam", package="Rsamtools")
to <- tempfile()
zipped <- bgzip(from, to)
idx <- indexTabix(zipped, "sam")

tab <- TabixFile(zipped, idx) 
}

\keyword{ manip }
