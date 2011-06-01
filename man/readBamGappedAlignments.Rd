\name{readBamGappedAlignments}
\Rdversion{1.1}

\alias{readBamGappedAlignments}
\alias{readBamGappedAlignments,character-method}
\alias{readBamGappedReads}
\alias{readBamGappedReads,character-method}

\title{Reading GappedAlignments or GappedReads objects from a BAM file}

\description{
Read a BAM file as a \link[GenomicRanges]{GappedAlignments}
or \link[ShortRead]{GappedReads} object.
}

\usage{
readBamGappedAlignments(file, index=file, use.names=FALSE, ..., which)
readBamGappedReads(file, index=file, use.names=FALSE, ..., which)
}

\arguments{

  \item{file}{The character(1) file name of the `BAM' file to be processed.}
  
  \item{index}{The character(1) name of the index file of the 'BAM' file
	being processed; this is given \emph{without} the '.bai' extension.}

  \item{use.names}{Construct the names of the returned object
     from the query template names (QNAME field)? If not (the default),
     then the returned object has no names.}

  \item{\dots}{Additional arguments, currently unused.}
  
  \item{which}{An object passed to \code{which} in
    \code{\linkS4class{ScanBamParam}} to specify ranges from which
    alignments will be retrieved. Valid types are described on the
    \code{\linkS4class{ScanBamParam}} help page.}

}

\details{
  
  See \code{?\link[GenomicRanges]{GappedAlignments-class}} for a
  description of \link[GenomicRanges]{GappedAlignments} objects.

  See \code{?\link[ShortRead]{GappedReads-class}} for a
  description of \link[ShortRead]{GappedReads} objects.

  See \code{?\link{scanBam}} for a description of the arguments. Unlike
  SAM/BAM records, we don't support unaligned queries so we discard
  those records.
  
}

\author{H. Pages}

\seealso{
  \link[GenomicRanges]{GappedAlignments-class},
  \link[ShortRead]{GappedReads-class},
  \code{\link{scanBam}}
}

\examples{
ex1_bamfile <- system.file("extdata", "ex1.bam", package="Rsamtools")
galn1 <- readBamGappedAlignments(ex1_bamfile)
galn1
names(galn1)
galn2 <- readBamGappedAlignments(ex1_bamfile, use.names=TRUE)
galn2
head(names(galn2))

greads1 <- readBamGappedReads(ex1_bamfile)
greads1
}

\keyword{methods}
\keyword{classes}
