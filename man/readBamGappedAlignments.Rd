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
readBamGappedAlignments(file, index=file, use.names=FALSE, param=NULL)
readBamGappedReads(file, index=file, use.names=FALSE, param=NULL)
}

\arguments{

  \item{file}{The character(1) file name of the `BAM' file to be processed.}
  
  \item{index}{The character(1) name of the index file of the 'BAM' file
    being processed; this is given \emph{without} the '.bai' extension.}

  \item{use.names}{Use the query template names (QNAME field) as the
    names of the returned object? If not (the default), then the returned
    object has no names.}

  \item{param}{\code{NULL} or an instance of \code{\linkS4class{ScanBamParam}}.
    Like for \code{\link{scanBam}}, this influences what fields and which
    records are imported. However, please note that the fields specified
    thru this \code{\linkS4class{ScanBamParam}} object will be loaded
    \emph{in addition} to any field required for generating the returned
    object (\link[GenomicRanges]{GappedAlignments} or
    \link[GenomicRanges]{GappedReads} object) and will be stored in its
    elementMetadata part.
    By default (i.e. \code{param=NULL}), no additional field is loaded and
    the flag used is \code{scanBamFlag(isUnmappedQuery=FALSE,
    isDuplicate=FALSE)} (i.e. all the records corresponding to mapped reads
    that are not PCR or optical duplicates are loaded).}
}

\details{
  
  See \code{?\link[GenomicRanges]{GappedAlignments-class}} for a
  description of \link[GenomicRanges]{GappedAlignments} objects.

  See \code{?\link[ShortRead]{GappedReads-class}} for a
  description of \link[ShortRead]{GappedReads} objects.

  See \code{?\link{scanBam}} for a description of the arguments.

}

\value{
  A \link[GenomicRanges]{GappedAlignments} object for
  \code{readBamGappedAlignments}.

  A \link[ShortRead]{GappedReads} object for \code{readBamGappedReads}.
}

\note{
  BAM records corresponding to unmapped reads or to reads that are PCR
  or optical duplicates are always ignored.
}

\author{H. Pages}

\seealso{
  \link[GenomicRanges]{GappedAlignments-class},
  \link[ShortRead]{GappedReads-class},
  \code{\link{scanBam}},
  \code{\link{ScanBamParam}}
}

\examples{
## Simple use:
bamfile <- system.file("extdata", "ex1.bam", package="Rsamtools")
galn1 <- readBamGappedAlignments(bamfile)
galn1
names(galn1)

## Using the 'use.names' arg:
galn2 <- readBamGappedAlignments(bamfile, use.names=TRUE)
galn2
head(names(galn2))

## Using the 'param' arg to load additional BAM fields:
param <- ScanBamParam(what=c("qual", "flag"))
galn3 <- readBamGappedAlignments(bamfile, param=param)
galn3
elementMetadata(galn3)

## Using the 'param' arg to load reads from particular regions.
## Note that if we don't specify 'what=character(0)' here, all
## the BAM fields will be loaded:
which <- RangesList(seq1=IRanges(1000, 2000),
                    seq2=IRanges(c(100, 1000), c(1000, 2000)))
param <- ScanBamParam(what=character(0), which=which)
galn4 <- readBamGappedAlignments(bamfile, param=param)
galn4

## Note that a given record is loaded one time for each region it
## belongs to (this is a scanBam() feature, readBamGappedAlignments()
## is based on scanBam()):
which <- IRangesList(seq2=IRanges(c(1563, 1567), width=1))
param <- ScanBamParam(which=which)
galn5 <- readBamGappedAlignments(bamfile, param=param)
galn5

## readBamGappedReads():
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
