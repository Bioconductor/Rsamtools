\name{seqnamesTabix}
\Rdversion{1.1}

\alias{seqnamesTabix}
\alias{seqnamesTabix,character-method}

\title{

  Retrieve sequence names defined in a tabix file.

}
\description{

  This function queries a tabix file, returning the names of the
  \sQuote{sequences} used as a key when creating the file.

}
\usage{

seqnamesTabix(file, ...)
\S4method{seqnamesTabix}{character}(file, ...)

}

\arguments{

  \item{file}{A \code{character(1)} file path or
    \code{\linkS4class{TabixFile}} instance pointing to a \sQuote{tabix}
    file.}

  \item{...}{Additional arguments, currently ignored.}

}

\value{

  A \code{character()} vector of sequence names present in the file.

}

\author{

  Martin Morgan <mtmorgan@fhcrc.org>.

}

\examples{
fl <- system.file("extdata", "example.gtf.gz", package="Rsamtools")
seqnamesTabix(fl)
}

\keyword{ manip }
