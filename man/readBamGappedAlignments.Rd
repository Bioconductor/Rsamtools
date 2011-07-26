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

  \item{use.names}{Use the query template names (QNAME field) as the
    names of the returned object? If not (the default), then the returned
    object has no names.}

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

  See \code{?\link{scanBam}} for a description of the arguments.

}

\note{
  BAM records corresponding to unaligned queries are ignored.
}

\author{H. Pages}

\seealso{
  \link[GenomicRanges]{GappedAlignments-class},
  \link[ShortRead]{GappedReads-class},
  \code{\link{scanBam}}
}

\examples{
bamfile <- system.file("extdata", "ex1.bam", package="Rsamtools")
galn1 <- readBamGappedAlignments(bamfile)
galn1
names(galn1)
galn2 <- readBamGappedAlignments(bamfile, use.names=TRUE)
galn2
head(names(galn2))

greads1 <- readBamGappedReads(bamfile)
greads1
names(greads1)
qseq(greads1)
greads2 <- readBamGappedReads(bamfile, use.names=TRUE)
greads2
head(names(greads2))
}

\keyword{methods}
\keyword{classes}
