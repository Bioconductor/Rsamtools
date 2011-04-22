\name{headerTabix}
\Rdversion{1.1}

\alias{headerTabix}
\alias{headerTabix,character-method}

\title{

  Retrieve sequence names defined in a tabix file.

}
\description{

  This function queries a tabix file, returning the names of the
  \sQuote{sequences} used as a key when creating the file.

}
\usage{

headerTabix(file, ...)
\S4method{headerTabix}{character}(file, ...)

}

\arguments{

  \item{file}{A \code{character(1)} file path or
    \code{\linkS4class{TabixFile}} instance pointing to a \sQuote{tabix}
    file.}

  \item{...}{Additional arguments, currently ignored.}

}

\value{

  A \code{list(4)} of the sequence names, column indicies used to sort
  the file, the number of lines skipped while indexing, and the comment
  character used while indexing.

}

\author{

  Martin Morgan <mtmorgan@fhcrc.org>.

}

\examples{
fl <- system.file("extdata", "example.gtf.gz", package="Rsamtools")
headerTabix(fl)
}

\keyword{ manip }
