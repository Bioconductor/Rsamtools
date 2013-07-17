\name{stackStringsFromBam}

\alias{stackStringsFromBam}


\title{Stack the read sequences stored in a BAM file on a region of interest}

\description{
  \code{stackStringsFromBam} lays the read sequences (or their quality
  strings) stored in a BAM file alongside the reference space, and stacks
  them on the specified region.
}

\usage{
stackStringsFromBam(file, param, use.names=FALSE, what="seq",
                    D.letter="-", N.letter="-",
                    Lpadding.letter="+", Rpadding.letter="+")
}

\arguments{
  \item{file}{
  }
  \item{param}{
  }
  \item{use.names}{
  }
  \item{what}{
  }
  \item{D.letter, N.letter}{
  }
  \item{Lpadding.letter, Rpadding.letter}{
  }
}

\details{
}

\value{
  An \link{XStringSet} object of the same class and length as \code{x}.
}

\author{H. Pages}

\seealso{
  \itemize{
    \item The \code{\link{readGAlignmentsFromBam}} function for loading read
          sequences (or their quality strings) from a BAM file (via a
          \link[GenomicRanges]{GAlignments} object).

    \item The \code{\link{sequenceLayer}} function for laying read sequences
          alongside the reference space, using their CIGARs.

    \item The \code{\link[Biostrings]{stackStrings}} function in the
          \pkg{Biostrings} package for stacking an arbitrary
          \link[Biostrings]{XStringSet} object.
  }
}

\examples{
## ---------------------------------------------------------------------
## A. EXAMPLE WITH TOY DATA
## ---------------------------------------------------------------------

bamfile <- BamFile(system.file("extdata", "ex1.bam", package="Rsamtools"))

stackStringsFromBam(bamfile, GRanges("seq1", IRanges(1, 60)))

options(showHeadLines=18)
options(showTailLines=6)
stackStringsFromBam(bamfile, GRanges("seq1", IRanges(61, 120)))

stacked_reads <- stackStringsFromBam(bamfile, "seq2:1509-1519")
stacked_reads  # deletion in read 13

stackStringsFromBam(bamfile, "seq2:1509-1519", what="qual")
consensusMatrix(stacked_reads)

## ---------------------------------------------------------------------
## B. EXAMPLE WITH REAL DATA
## ---------------------------------------------------------------------

library(RNAseqData.HNRNPC.bam.chr14)
bamfile <- BamFile(RNAseqData.HNRNPC.bam.chr14_BAMFILES[1])

## My Region Of Interest:
my_ROI <- GRanges("chr14", IRanges(19650095, 19650159))

readGAlignments(bamfile, param=ScanBamParam(which=my_ROI))
stackStringsFromBam(bamfile, my_ROI)
}

\keyword{methods}
\keyword{manip}
