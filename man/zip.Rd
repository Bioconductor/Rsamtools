\name{Compression}
\Rdversion{1.1}

\alias{bgzip}
\alias{razip}

\title{

  File compression for tabix (bgzip) and fasta (razip) files.

}
\description{

  These functions compress files for use in other parts of
  \pkg{Rsamtools}: \code{bgzip} for tabix files, \code{razip} for
  random-access fasta files.

}
\usage{

bgzip(file, dest=sprintf("\%s.gz", file), overwrite = FALSE)
razip(file, dest=sprintf("\%s.rz", file), overwrite = FALSE)

}

\arguments{

  \item{file}{A character(1) path to an existing file. This file will be
    compressed.}

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

from <- system.file("extdata", "ex1.sam", package="Rsamtools")
to <- tempfile()
zipped <- bgzip(from, to)

}

\keyword{ manip }
