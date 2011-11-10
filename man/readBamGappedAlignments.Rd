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

  \item{file}{The path to the `BAM' file to read.}
  
  \item{index}{The path to the index file of the 'BAM' file
    to read; this is given \emph{without} the '.bai' extension.}

  \item{use.names}{Use the query template names (QNAME field) as the
    names of the returned object? If not (the default), then the returned
    object has no names.}

  \item{param}{\code{NULL} or an instance of \code{\linkS4class{ScanBamParam}}.
    Like for \code{\link{scanBam}}, this influences what fields and which
    records are imported. However, please note that the fields specified
    thru this \code{\linkS4class{ScanBamParam}} object will be loaded
    \emph{in addition} to any field required for generating the returned
    object (\link[GenomicRanges]{GappedAlignments} or
    \link[ShortRead]{GappedReads} object) and will be stored in its
    elementMetadata part.
    By default (i.e. \code{param=NULL}), no additional field is loaded and
    the flag used is \code{scanBamFlag(isUnmappedQuery=FALSE)} (i.e. records
    corresponding to unmapped reads are not loaded).}
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
  BAM records corresponding to unmapped reads are always ignored.

  Starting with Rsamtools 1.7.1 (BioC 2.10), PCR or optical duplicates
  are loaded by default (use \code{scanBamFlag(isDuplicate=FALSE)} to
  drop them).
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

## Using the 'param' arg to drop PCR or optical duplicates and load
## additional BAM fields:
param <- ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE),
                      what=c("qual", "flag"))
galn3 <- readBamGappedAlignments(bamfile, param=param)
galn3
elementMetadata(galn3)

## Using the 'param' arg to load reads from particular regions.
## Note that if we weren't providing a 'what' argument here, all the
## BAM fields would be loaded:
which <- RangesList(seq1=IRanges(1000, 2000),
                    seq2=IRanges(c(100, 1000), c(1000, 2000)))
param <- ScanBamParam(which=which)
galn4 <- readBamGappedAlignments(bamfile, param=param)
galn4

## Note that a given record is loaded one time for each region it
## belongs to (this is a scanBam() feature, readBamGappedAlignments()
## is based on scanBam()):
which <- IRangesList(seq2=IRanges(c(1563, 1567), width=1))
param <- ScanBamParam(which=which)
galn5 <- readBamGappedAlignments(bamfile, param=param)
galn5

## Using the 'param' arg to load tags. Except for MF and Aq, the tags
## specified below are predefined tags (see the SAM Spec for the list
## of predefined tags and their meaning).
param <- ScanBamParam(tag=c("MF", "Aq", "NM", "UQ", "H0", "H1"),
                      what="isize")
galn6 <- readBamGappedAlignments(bamfile, param=param)
elementMetadata(galn6)  # "tag" cols always after "what" cols

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
