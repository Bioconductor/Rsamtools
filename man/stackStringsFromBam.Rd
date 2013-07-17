\name{stackStringsFromBam}

\alias{stackStringsFromBam}


\title{Stack the read sequences stored in a BAM file on a region of interest}

\description{
  \code{stackStringsFromBam} lays the read sequences (or their quality
  strings) stored in a BAM file alongside the reference space, and stacks
  them on the specified region.
}

\usage{
stackStringsFromBam(file, index=file, param,
                    what="seq", use.names=FALSE,
                    D.letter="-", N.letter="-",
                    Lpadding.letter="+", Rpadding.letter="+")
}

\arguments{
  \item{file, index}{
    The path to the BAM file to read, and to the index file of the BAM file
    to read, respectively. The latter is given \emph{without} the '.bai'
    extension. See \code{\link{scanBam}} for more information.
  }
  \item{param}{
    A \link{ScanBamParam} object containing exactly 1 genomic region
    (i.e. \code{unlist(bamWhich(param))} must have length 1).

    Alternatively, \code{param} can be a \link[GenomicRanges]{GRanges} or
    \link[IRanges]{RangesList} object containing exactly 1 genomic region,
    or a character string specifying a single genomic region (in the
    \code{"chr14:5201-5300"} format).
  }
  \item{what}{
    A single string. Either \code{"seq"} or \code{"qual"}. If \code{"seq"},
    the read sequences will be stacked. If \code{"qual"}, the read quality
    strings will be stacked.
  }
  \item{use.names}{
    Use the query template names (QNAME field) as the names of the returned
    object? If not (the default), then the returned object has no names.
  }
  \item{D.letter, N.letter}{
    A single letter used as a filler for injections. The 2 arguments are
    passed down to the \code{\link{sequenceLayer}} function.
    See \code{?\link{sequenceLayer}} for more details.
  }
  \item{Lpadding.letter, Rpadding.letter}{
    A single letter to use for padding the sequences on the left, and another
    one to use for padding on the right. The 2 arguments are passed down to
    the \code{\link[Biostrings]{stackStrings}} function defined in the
    \pkg{Biostrings} package.
    See \code{?\link[Biostrings]{stackStrings}} in the \pkg{Biostrings}
    package for more details.
  }
}

\details{
  \code{stackStringsFromBam} performs the 3 following steps:
  \enumerate{
    \item Load the read sequences (or their quality strings) from the BAM
          file. Only the read sequences that overlap with the specified region
          are loaded. This is done with the
          \code{\link{readGAlignmentsFromBam}} function. Note that if the
          file contains paired-end reads, the pairing is ignored.

    \item Lay the sequences alongside the reference space, using their CIGARs.
          This is done with the \code{\link{sequenceLayer}} function.

    \item Stack them on the specified region. This is done with the
          \code{\link[Biostrings]{stackStrings}} function defined in the
          \pkg{Biostrings} package.
  }
}

\note{
  TWO IMPORTANT CAVEATS:

  Specifying a big genomic region, say >= 100000 bp, can require a lot of
  memory (especially with high coverage reads) and is not recommended.

  Paired-end reads are treated as single-end reads (i.e. they're not paired).
}

\value{
  A rectangular (i.e. constant-width) \link[Biostrings]{DNAStringSet} object
  if \code{what} is \code{"seq"} or \link[Biostrings]{BStringSet} object
  if \code{what} is \code{"qual"}.
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

stackStringsFromBam(bamfile, param=GRanges("seq1", IRanges(1, 60)))

options(showHeadLines=18)
options(showTailLines=6)
stackStringsFromBam(bamfile, param=GRanges("seq1", IRanges(61, 120)))

stacked_qseq <- stackStringsFromBam(bamfile, param="seq2:1509-1519")
stacked_qseq  # deletion in read 13

stackStringsFromBam(bamfile, param="seq2:1509-1519", what="qual")
consensusMatrix(stacked_qseq)

## ---------------------------------------------------------------------
## B. EXAMPLE WITH REAL DATA
## ---------------------------------------------------------------------

library(RNAseqData.HNRNPC.bam.chr14)
bamfile <- BamFile(RNAseqData.HNRNPC.bam.chr14_BAMFILES[1])

## My Region Of Interest:
my_ROI <- GRanges("chr14", IRanges(19650095, 19650159))

readGAlignments(bamfile, param=ScanBamParam(which=my_ROI))
stackStringsFromBam(bamfile, param=my_ROI)
}

\keyword{methods}
\keyword{manip}
