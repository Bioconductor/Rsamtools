\name{cigar-utils}
\Rdversion{1.1}
\alias{cigar-utils}
\alias{cigarToReadWidth}
\alias{cigarToIRanges}
\alias{cigarToIRangesList}
\alias{cigarToCigarTable}
\alias{summarizeCigarTable}

\title{
  CIGAR utility functions
}
\description{
  Utility functions for low-level CIGAR manipulation.
}
\usage{
cigarToReadWidth(cigar, after.hard.clipping=FALSE)

cigarToIRanges(cigar, drop.D.ranges=FALSE, merge.ranges=TRUE)

cigarToIRangesList(cigar, rname, pos, flag=NA,
                   drop.D.ranges=FALSE, merge.ranges=TRUE)

cigarToCigarTable(cigar)

summarizeCigarTable(x)
}

\arguments{
  \item{cigar}{
    A character vector/factor containing the extended CIGAR
    string for each read.
  }
  \item{after.hard.clipping}{
    Should the returned widths be those of the reads before
    or after "hard clipping"? Hard clipping of a read is
    encoded with an H in the CIGAR.
    If NO (\code{after.hard.clipping=FALSE}, the default),
    then the returned widths are those of the original reads.
    If YES (\code{after.hard.clipping=TRUE}), then the
    returned widths are the lengths of the query sequences
    stored in the SAM/BAM file.
  }
  \item{drop.D.ranges}{
    Should the ranges corresponding to a deletion from the
    reference (encoded with a D in the CIGAR) be dropped?
    By default we keep them to be consistent with the pileup tool
    from SAMtools.
    Note that, when \code{drop.D.ranges} is \code{TRUE}, then Ds
    and Ns in the CIGAR are equivalent.
  }
  \item{merge.ranges}{
    Should adjacent ranges coming from the same cigar be merged or not?
    Using \code{TRUE} (the default) can significantly reduce the
    size of the returned object.
  }
  \item{rname}{
    A character vector/factor containing the name of the
    reference sequence associated with each read (i.e. the
    name of the sequence the read has been aligned to).
  }
  \item{pos}{
    An integer vector containing the 1-based leftmost
    position/coordinate for each (eventually clipped) read
    sequence.
  }
  \item{flag}{
    \code{NA} or an integer vector containing the SAM flag for
    each read.
    According to the SAM specs, when bit 0x400 in the flag
    is ON then "the read is either a PCR duplicate or an optical
    duplicate". Therefore when \code{flag} is provided, reads
    corresponding to a flag where this bit is ON are ignored.
  }
  \item{x}{
    A \link[IRanges]{DataFrame} produced by \code{cigarToCigarTable}.
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

  For \code{cigarToCigarTable}: A frequency table of the CIGARs
  in the form of a \link[IRanges]{DataFrame} with two columns:
  \code{cigar} (a \link[IRanges]{CompressedRleList}) and \code{count}
  (an integer).

  For \code{summarizeCigarTable}: A list with two elements:
  \code{AlignedCharacters} (integer) and \code{Indels} (matrix)
}

\references{
  \url{http://samtools.sourceforge.net/}
}

\author{
  H. Pages
}

\seealso{
  \link[IRanges]{IRanges-class},
  \link[IRanges]{IRangesList-class},
  \code{\link[IRanges]{coverage}},
  \link[IRanges]{RleList-class}
}

\examples{
## ---------------------------------------------------------------------
## A. SIMPLE EXAMPLES
## ---------------------------------------------------------------------

cigarToReadWidth("3H15M55N4M2I6M2D5M6S")
cigarToReadWidth("3H15M55N4M2I6M2D5M6S", after.hard.clipping=TRUE)

cigarToIRanges("3H15M55N4M2I6M2D5M6S")
cigarToIRanges("3H15M55N4M2I6M2D5M6S", merge.ranges=FALSE)
cigarToIRanges("3H15M55N4M2I6M2D5M6S", drop.D.ranges=TRUE)

cigar <- c("40M", "3H15M55N4M2I6M2D5M6S", "2S10M2000N15M", "3H25M5H")
rname <- c("chr6", "chr6", "chr2", "chr6")
pos <- c(1, 1001, 1,  351)
cigarToIRangesList(cigar, rname, pos)

cigarToCigarTable(cigar)
cigarToCigarTable(cigar)[,"cigar"]
cigarToCigarTable(cigar)[,"count"]

summarizeCigarTable(cigarToCigarTable(cigar))

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
  system.time(irl <- cigarToIRangesList(cigar, rname, pos))

  ## Internally, cigarToIRangesList() turns 'rname' into a factor
  ## before starting the calculation. Hence it will run sligthly
  ## faster if 'rname' is already a factor.
  rname2 <- as.factor(rname)
  system.time(irl2 <- cigarToIRangesList(cigar, rname2, pos))

  ## The size of the resulting IRangesList object is about 160M:
  object.size(irl)
}

## ---------------------------------------------------------------------
## C. COMPUTE THE COVERAGE OF THE READS STORED IN A BAM FILE
## ---------------------------------------------------------------------
## The information stored in a BAM file can be used to compute the
## "coverage" of the mapped reads i.e. the number of reads that hit any
## given position in the reference genome.
## The following function takes the path to a BAM file and returns an
## object representing the coverage of the mapped reads stored in the
## file. The returned object is an RleList object named with the names
## of the reference sequences that actually receive some coverage.

extractCoverageFromBAM <- function(file)
{
  ## This ScanBamParam object allows us to load only the necessary
  ## information from the file.
  param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE,
                                         isDuplicate=FALSE),
                        what=c("rname", "pos", "cigar"))
  bam <- scanBam(file, param=param)[[1]]
  irl <- cigarToIRangesList(as.vector(bam$cigar@.cigar), bam$rname, bam$pos)
  irl <- irl[elementLengths(irl) != 0]  # drop empty elements
  coverage(irl)
}


fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
extractCoverageFromBAM(fl)
}

\keyword{manip}
