\name{Rsamtools-package}
\Rdversion{1.1}
\alias{Rsamtools-package}
\alias{Rsamtools}
\alias{bzfile-class}
\alias{fifo-class}
\alias{gzfile-class}
\alias{pipe-class}
\alias{unz-class}
\alias{url-class}

\docType{package}
\title{
  'samtools' aligned sequence utilities interface
}
\description{
  This package provides facilities for parsing samtools BAM (binary)
  files representing aligned sequences.
}
\details{

  See \code{packageDescription('Rsamtools')} for package details. A
  useful starting point is the \code{\link{scanBam}} manual page.

}
\author{
Author: Martin Morgan

Maintainer: Biocore Team c/o BioC user list <bioconductor@stat.math.ethz.ch>
}
\references{

  The current source code for samtools and bcftools is from
  \url{https://github.com/samtools/samtools}.  Additional material is at
  \url{http://samtools.sourceforge.net/}.

}
\note{

  This package documents the following classes for purely internal
  reasons, see help pages in other packages: \code{bzfile}, \code{fifo},
  \code{gzfile}, \code{pipe}, \code{unz}, \code{url}.

}
\keyword{ package }
\examples{
packageDescription('Rsamtools')
}
