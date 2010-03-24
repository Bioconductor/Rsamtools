\name{readBamGappedAlignments}
\Rdversion{1.1}

\alias{readBamGappedAlignments}
\alias{readBamGappedAlignments,character-method}

\title{GappedAlignments objects}

\description{
Read a BAM file as a \link[GenomicRanges]{GappedAlignments} object.
}

\usage{
readBamGappedAlignments(file, index, ..., which, ans.subtype="Alignments0")
}

\arguments{

  \item{file}{The character(1) file name of the `BAM' file to be processed.}
  
  \item{index}{The character(1) name of the index file of the 'BAM' file
	being processed; this is given \emph{without} the '.bai' extension.}

  \item{\dots}{Additional arguments, currently unused.}
  
  \item{which}{An object passed to \code{which} in
    \code{\linkS4class{ScanBamParam}} to specify ranges from which
    alignments will be retrieved. Valid types are described on the
    \code{\linkS4class{ScanBamParam}} help page.}

  \item{ans.subtype}{A (temporary) flag indicating the return type of
    the function.}

}

\details{
  
  See \code{?\link[GenomicRanges]{GappedAlignments-class}} for a
  description of \link[GenomicRanges]{GappedAlignments} objects.

  See \code{?\link{scanBam}} for a description of the arguments. Unlike
  SAM/BAM records, we don't support unaligned queries so we discard
  those records.
  
}

\author{H. Pages}

\seealso{
  \link[GenomicRanges]{GappedAlignments-class},
  \code{\link{scanBam}}
}

\examples{
aln1_file <- system.file("extdata", "ex1.bam", package="Rsamtools")
aln1 <- readBamGappedAlignments(aln1_file)
aln1
}

\keyword{methods}
\keyword{classes}
