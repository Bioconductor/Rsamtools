\name{cigar-utils}
\Rdversion{1.1}
\alias{cigar-utils}
\alias{cigarToIRangesList}

\title{
  CIGAR utility functions
}
\description{
  Utility functions for low-level CIGAR manipulation.
}
\usage{
cigarToIRangesList(rname, strand, pos, cigar)
}

\arguments{
  \item{rname}{
    A character vector containing the name of the reference sequence
    associated with each read (i.e. the name of the sequence the read
    has been aligned to).
  }
  \item{strand}{
    Ignored for now.
  }
  \item{pos}{
    An integer vector containing the 1-based leftmost position/coordinate
    of the clipped read sequence.
  }
  \item{cigar}{
    A character vector containing the extended CIGAR string for each
    read.
  }
}

\value{
  An \link[IRanges]{IRangesList} object.
}

\references{
  \url{http://samtools.sourceforge.net/}
}

\author{
  H. Pages
}

\examples{
rname <- c("chr6", "chr6", "chr2", "chr6")
pos <- c(1, 300, 21,  320)
cigar <- c("9M100D15M", "99M", "999M4M88I5M3D15M", "5M")
cigarToIRangesList(rname, NULL, pos, cigar)
}

\keyword{ manip }
