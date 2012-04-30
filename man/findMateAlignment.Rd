\name{findMateAlignment}

\alias{findMateAlignment}
\alias{findMateAlignment2}
\alias{makeGappedAlignmentPairs}

\title{Pairing the elements of a GappedAlignments object}

\description{
Utilities for pairing the elements of a \link[GenomicRanges]{GappedAlignments}
object.
}

\usage{
findMateAlignment(x, verbose=FALSE)
makeGappedAlignmentPairs(x, use.names=FALSE, keep.cols=NULL)
}

\arguments{
  \item{x}{
    A \link[GenomicRanges]{GappedAlignments} object.
  }
  \item{verbose}{
    If \code{TRUE}, then \code{findMateAlignment} will print some
    details about what is currently going on. Mostly useful for debugging.
  }
  \item{use.names}{
    Whether the names on the input object should be propagated to the
    returned object or not.
  }
  \item{keep.cols}{
    Names of the metadata columns to keep in the returned
    \link{GappedAlignmentPairs} object.
  }
}

\details{
  COMING SOON...
}

\value{
  COMING SOON...
}

\author{H. Pages}

\seealso{
  \link[GenomicRanges]{GappedAlignments-class},
  \link[GenomicRanges]{GappedAlignmentPairs-class},
  \code{\link{readBamGappedAlignments}},
  \code{\link{readBamGappedAlignmentPairs}}
}

\examples{
## COMING SOON...
}

\keyword{methods}
\keyword{classes}
