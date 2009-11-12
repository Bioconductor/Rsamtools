\name{cigar-utils}
\Rdversion{1.1}
\alias{cigar-utils}
\alias{cigarToReadWidth}
\alias{cigarToIRanges}
\alias{cigarToIRangesList}

\title{
  CIGAR utility functions
}
\description{
  Utility functions for low-level CIGAR manipulation.
}
\usage{
cigarToReadWidth(cigar)
cigarToIRanges(cigar)
cigarToIRangesList(cigar, rname, strand, pos)
}

\arguments{
  \item{cigar}{
    A character vector/factor containing the extended CIGAR
    string for each read.
  }
  \item{rname}{
    A character vector/factor containing the name of the
    reference sequence associated with each read (i.e. the
    name of the sequence the read has been aligned to).
  }
  \item{strand}{
    Ignored for now.
  }
  \item{pos}{
    An integer vector containing the 1-based leftmost
    position/coordinate for each (eventually clipped) read
    sequence.
  }
}

\value{
  For \code{cigarToReadWidth}: An integer vector of the same length
  as \code{cigar} where each element is the width of the read as
  inferred from the corresponding element in \code{cigar}.

  For \code{cigarToIRanges}: An \link[IRanges]{IRanges} object
  describing where the bases in the read align with respect to
  an imaginary reference sequence assuming that the leftmost
  aligned base is at position 1 in the reference (i.e. at the
  first position).

  For \code{cigarToIRangesList}: A named \link[IRanges]{IRangesList}
  object with one element (\link[IRanges]{IRanges}) per unique reference
  sequence.
}

\references{
  \url{http://samtools.sourceforge.net/}
}

\author{
  H. Pages
}

\examples{
## ---------------------------------------------------------------------
## A. SIMPLE EXAMPLES
## ---------------------------------------------------------------------

cigarToReadWidth("15M2I6M2D5M6S")

cigarToIRanges("15M2I6M2D5M6S")

cigar <- c("40M", "15M2I6M2D5M6S", "10M2000N15M", "3H25M5H")
rname <- c("chr6", "chr6", "chr2", "chr6")
pos <- c(1, 1001, 1,  351)
cigarToIRangesList(cigar, rname, NULL, pos)

## ---------------------------------------------------------------------
## B. PERFORMANCE
## ---------------------------------------------------------------------

if (interactive()) {
  ## We simulate 20 millions aligned reads, all 40-mers. 95% of them
  ## align with no indels. 5% align with a big deletion in the
  ## reference. In the context of an RNAseq experiment, those 5% would
  ## be suspected to be "junction reads".
  set.seed(123)
  nreads <- 20000000L
  njunctionreads <- nreads * 5L / 100L
  cigar <- character(nreads)
  cigar[] <- "40M"
  junctioncigars <- paste(
      paste(10:30, "M", sep=""),
      paste(sample(80:8000, njunctionreads, replace=TRUE), "N", sep=""),
      paste(30:10, "M", sep=""), sep="")
  cigar[sample(nreads, njunctionreads)] <- junctioncigars
  some_fake_rnames <- paste("chr", c(1:6, "X"), sep="")
  rname <- sample(some_fake_rnames, nreads, replace=TRUE)
  pos <- sample(80000000L, nreads, replace=TRUE)

  ## The following takes < 10 sec. to complete:
  system.time(irl <- cigarToIRangesList(cigar, rname, NULL, pos))

  ## Internally, cigarToIRangesList() turns 'rname' into a factor
  ## before starting the calculation. Hence it will run sligthly
  ## faster if 'rname' is already a factor.
  rname2 <- as.factor(rname)
  system.time(irl2 <- cigarToIRangesList(cigar, rname2, NULL, pos))

  ## The size of the resulting IRangesList object is about 160M:
  object.size(irl)
}
}

\keyword{manip}
