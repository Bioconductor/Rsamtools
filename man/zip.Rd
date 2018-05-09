\name{Compression}
\Rdversion{1.1}

\alias{bgzip}
\alias{razip}

\title{

  File compression for tabix (bgzip) and fasta (razip) files.

  IMPORTANT NOTE: Starting with Rsamtools 1.99.0 (Bioconductor 3.9),
  \code{razip()} is defunct. Please use \code{bgzip()} instead.
}
\description{

  \code{bgzip} compresses tabix (e.g. SAM or VCF) or FASTA files
  to a format that allows indexing for later fast random-access.

}
\usage{

bgzip(file, dest=sprintf("\%s.bgz", sub("\\\\.gz$", "", file)),
      overwrite = FALSE)

## Defunct!
razip(file, dest=sprintf("\%s.rz", sub("\\\\.gz$", "", file)),
      overwrite = FALSE)

}

\arguments{

  \item{file}{A character(1) path to an existing uncompressed or
    gz-compressed file. This file will be compressed.}

  \item{dest}{A character(1) path to a file. This will be the compressed
    file. If \code{dest} exists, then it is only over-written when
    \code{overwrite=TRUE}.}

  \item{overwrite}{A logical(1) indicating whether \code{dest} should
    be over-written, if it already exists.}

}

\value{

  The full path to \code{dest}.

}
\references{
  \url{http://samtools.sourceforge.net/}
}
\author{

  Martin Morgan <mtmorgan@fhcrc.org>

}

\seealso{

  \code{\link{TabixFile}}, \code{\link{FaFile}}.

}

\examples{

from <- system.file("extdata", "ex1.sam", package="Rsamtools",
                    mustWork=TRUE)
to <- tempfile()
zipped <- bgzip(from, to)

}

\keyword{ manip }
