\name{TabixInput}
\Rdversion{1.1}

\alias{scanTabix}
\alias{scanTabix,character,GRanges-method}
\alias{scanTabix,character,RangedData-method}
\alias{scanTabix,character,RangesList-method}

\title{

  Operations on `tabix' (indexed, tab-delimited) files.

}
\description{
  Scan compressed, sorted, tabix-indexed, tab-delimited files.
}
\usage{

scanTabix(file, ..., param)
\S4method{scanTabix}{character,RangesList}(file, ..., param)
\S4method{scanTabix}{character,RangedData}(file, ..., param)
\S4method{scanTabix}{character,GRanges}(file, ..., param)

}

\arguments{

  \item{file}{The character() file name(s) of the tabix file be
    processed, or more flexibly an instance of class
    \code{\link{TabixFile}}.}

  \item{param}{A instance of \code{GRanges}, \code{RangedData},  or
    \code{RangesList} provide the sequence names and regions to be
    parsed.}


  \item{...}{Additional arguments, currently ignored.}

}

\value{

  \code{scanTabix} returns a list, with one element per region. Each element
  of the list is a character vector representing records in the region.

}

\references{

  \url{http://samtools.sourceforge.net/tabix.shtml}

}

\author{

  Martin Morgan <mtmorgan@fhcrc.org>.

}

\examples{
example(TabixFile)
}

\keyword{ manip }
