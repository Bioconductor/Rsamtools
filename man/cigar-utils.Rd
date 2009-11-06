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
    A character vector/factor containing the name of the
    reference sequence associated with each read (i.e. the
    name of the sequence the read has been aligned to).
  }
  \item{strand}{
    Ignored for now.
  }
  \item{pos}{
    An integer vector containing the 1-based leftmost
    position/coordinate of the clipped read sequence.
  }
  \item{cigar}{
    A character vector/factor containing the extended CIGAR
    string for each read.
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
## ---------------------------------------------------------------------
## A. SMALL EXAMPLES
## ---------------------------------------------------------------------
rname <- c("chr6", "chr6", "chr2", "chr6")
pos <- c(1, 300, 21,  320)
cigar <- c("9M100D15M", "99M", "999M4M88I5M3D15M", "5M")
cigarToIRangesList(rname, NULL, pos, cigar)

## ---------------------------------------------------------------------
## B. PERFORMANCE
## ---------------------------------------------------------------------
if (interactive()) {
  ## We simulate 20 millions aligned reads, all 40-mers. 95% of them
  ## align with no indels. 5% align with a big deletion (from the
  ## reference to the read). In the context of an RNAseq experiment,
  ## those 5% would be suspected to be "junction reads".
  nreads <- 20000000L
  njunctionreads <- nreads * 5L / 100L
  some_fake_rnames <- paste("chr", c(1:6, "X"), sep="")
  rname <- sample(some_fake_rnames, nreads, replace=TRUE)
  pos <- sample(80000000L, nreads, replace=TRUE)
  cigar <- character(nreads)
  cigar[] <- "40M"
  junctioncigars <- paste(
      paste(10:30, "M", sep=""),
      paste(sample(80:8000, njunctionreads, replace=TRUE), "N", sep=""),
      paste(30:10, "M", sep=""), sep="")
  cigar[sample(nreads, njunctionreads)] <- junctioncigars

  ## The following takes between 13 and 15 sec. to complete:
  system.time(irl <- cigarToIRangesList(rname, NULL, pos, cigar))

  ## Internally, cigarToIRangesList() turns 'rname' into a factor
  ## before starting the calculation. Hence it will run sligthly
  ## faster if 'rname' is already a factor.
  rname2 <- as.factor(rname)
  system.time(irl2 <- cigarToIRangesList(rname2, NULL, pos, cigar))

  ## The size of the resulting IRangesList object is about 160M:
  object.size(irl)
}
}

\keyword{manip}
