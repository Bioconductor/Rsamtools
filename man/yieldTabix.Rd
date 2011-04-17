\name{yieldTabix}
\Rdversion{1.1}

\alias{yieldTabix}

\title{

  Yield records from a stream of records contained in a file or other
  object.

}
\description{
  Return the next set of records from an already opened file.
}
\usage{

yieldTabix(file, ..., yieldSize=1000000L)

}

\arguments{

  \item{file}{An opened instance of a type for which a \code{yield}
    method exists.}

  \item{yieldSize}{The number of records to return on each invocation.}


  \item{...}{Additional arguments, currently ignored.}

}

\value{

  \code{yield} returns the next \code{yieldSize} records from
  \code{file}, in a format defined by the method.

}

\author{

  Martin Morgan <mtmorgan@fhcrc.org>.

}

\examples{
showMethods(yieldTabix)
example(TabixFile)
}

\keyword{ manip }
