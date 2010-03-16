\name{readBAMasGappedAlignments}
\Rdversion{1.1}

\alias{readBAMasGappedAlignments}
\alias{readBAMasGappedAlignments,character-method}

\title{GappedAlignments objects}

\description{
Read a BAM file as a \link[GenomicRanges]{GappedAlignments} object.
}

\usage{
readBAMasGappedAlignments(file, index, ..., which, ans.subtype="Alignments0")
}

\details{
See \code{?\link[GenomicRanges]{GappedAlignments-class}} for a
description of \link[GenomicRanges]{GappedAlignments} objects.

See \code{?\link{scanBam}} for a description of the arguments.
Unlike SAM/BAM records, we don't support unaligned queries so
we discard those records.
}

\author{H. Pages}

\seealso{
  \link[GenomicRanges]{GappedAlignments-class},
  \code{\link{scanBam}}
}

\examples{
aln1_file <- system.file("extdata", "ex1.bam", package="Rsamtools")
aln1 <- readBAMasGappedAlignments(aln1_file)
aln1
}

\keyword{methods}
\keyword{classes}
