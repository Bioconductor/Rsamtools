\name{readBamGappedAlignments}
\Rdversion{1.1}

\alias{readBamGappedAlignments}
\alias{readBamGappedAlignments,character-method}
\alias{readBamGappedReads}
\alias{readBamGappedReads,character-method}
\alias{readBamGappedAlignmentPairs}
\alias{readBamGappedAlignmentPairs,character-method}

\title{Reading a BAM file into a GappedAlignments, GappedReads, or
       GappedAlignmentPairs object}

\description{
Read a BAM file into a \link[GenomicRanges]{GappedAlignments},
\link[ShortRead]{GappedReads}, or \link[GenomicRanges]{GappedAlignmentPairs}
object.
}

\usage{
readBamGappedAlignments(file, index=file, use.names=FALSE, param=NULL)
readBamGappedReads(file, index=file, use.names=FALSE, param=NULL)
readBamGappedAlignmentPairs(file, index=file, use.names=FALSE, param=NULL)
}

\arguments{
  \item{file, index}{The path to the BAM file to read, and to the index
    file of the BAM file to read, respectively. The latter is given
    \emph{without} the '.bai' extension. See \code{\link{scanBam}} for more
    information.}

  \item{use.names}{Use the query template names (QNAME field) as the
    names of the returned object? If not (the default), then the returned
    object has no names.}

  \item{param}{\code{NULL} or an instance of \code{\linkS4class{ScanBamParam}}.
    Like for \code{\link{scanBam}}, this influences what fields and which
    records are imported. However, note that the fields specified
    thru this \code{\linkS4class{ScanBamParam}} object will be loaded
    \emph{in addition} to any field required for generating the returned
    object (\link[GenomicRanges]{GappedAlignments},
    \link[ShortRead]{GappedReads}, or
    \link[GenomicRanges]{GappedAlignmentPairs} object), but only the fields
    requested by the user will actually be kept as metadata columns of the
    object.

    By default (i.e. \code{param=NULL}), no additional field is loaded and
    the flag used is \code{scanBamFlag(isUnmappedQuery=FALSE)} for
    \code{readBamGappedAlignments} and \code{readBamGappedReads}
    (i.e. only records corresponding to mapped reads are loaded),
    and \code{scanBamFlag(isUnmappedQuery=FALSE, isPaired=TRUE,
    hasUnmappedMate=FALSE)} for \code{readBamGappedAlignmentPairs}
    (i.e. only records corresponding to paired-end reads with both ends
    mapped are loaded).}
}

\details{
  See \code{?\link[GenomicRanges]{GappedAlignments-class}} for a
  description of \link[GenomicRanges]{GappedAlignments} objects.

  See \code{?\link[ShortRead]{GappedReads-class}} for a
  description of \link[ShortRead]{GappedReads} objects.

  \code{readBamGappedAlignmentPairs} proceeds in 2 steps:
  \enumerate{
    \item Load the BAM file into a \link[GenomicRanges]{GappedAlignments}
          object with \code{readBamGappedAlignments};
    \item Turn this \link[GenomicRanges]{GappedAlignments} object into a
          \link[GenomicRanges]{GappedAlignmentPairs} object by pairing
          its elements.
  }
  See \code{?\link[GenomicRanges]{GappedAlignmentPairs-class}} for a
  description of \link[GenomicRanges]{GappedAlignmentPairs} objects,
  and \code{?\link{makeGappedAlignmentPairs}} for the details of the
  pairing procedure.
}

\value{
  A \link[GenomicRanges]{GappedAlignments} object for
  \code{readBamGappedAlignments}.

  A \link[ShortRead]{GappedReads} object for \code{readBamGappedReads}.

  A \link[GenomicRanges]{GappedAlignmentPairs} object for
  \code{readBamGappedAlignmentPairs}.
  Note that a BAM (or SAM) file can in theory contain a mix of single-end
  and paired-end reads, but in practise it seems that single-end and
  paired-end are not mixed. In other words, the value of flag bit 0x1
  (\code{isPaired}) is the same for all the records in a file.
  So if \code{readBamGappedAlignmentPairs} returns a
  \link[GenomicRanges]{GappedAlignmentPairs} object of length zero,
  this almost certainly means that the BAM (or SAM) file contains
  alignments for single-end reads (although it could also mean that the
  user-supplied \code{\linkS4class{ScanBamParam}} is filtering out everything,
  or that the file is empty, or that all the records in the file correspond
  to unmapped reads).
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
  \link[GenomicRanges]{GappedAlignmentPairs-class},
  \code{\link{makeGappedAlignmentPairs}},
  \code{\link{scanBam}},
  \code{\link{ScanBamParam}}
}

\examples{
## ---------------------------------------------------------------------
## A. readBamGappedAlignments()
## ---------------------------------------------------------------------

## Simple use:
bamfile <- system.file("extdata", "ex1.bam", package="Rsamtools")
gal1 <- readBamGappedAlignments(bamfile)
gal1
names(gal1)

## Using the 'use.names' arg:
gal2 <- readBamGappedAlignments(bamfile, use.names=TRUE)
gal2
head(names(gal2))

## Using the 'param' arg to drop PCR or optical duplicates and load
## additional BAM fields:
param <- ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE),
                      what=c("qual", "flag"))
gal3 <- readBamGappedAlignments(bamfile, param=param)
gal3
mcols(gal3)

## Using the 'param' arg to load reads from particular regions.
## Note that if we weren't providing a 'what' argument here, all the
## BAM fields would be loaded:
which <- RangesList(seq1=IRanges(1000, 2000),
                    seq2=IRanges(c(100, 1000), c(1000, 2000)))
param <- ScanBamParam(which=which)
gal4 <- readBamGappedAlignments(bamfile, param=param)
gal4

## Note that a given record is loaded one time for each region it
## belongs to (this is a scanBam() feature, readBamGappedAlignments()
## is based on scanBam()):
which <- IRangesList(seq2=IRanges(c(1563, 1567), width=1))
param <- ScanBamParam(which=which)
gal5 <- readBamGappedAlignments(bamfile, param=param)
gal5

## Using the 'param' arg to load tags. Except for MF and Aq, the tags
## specified below are predefined tags (see the SAM Spec for the list
## of predefined tags and their meaning).
param <- ScanBamParam(tag=c("MF", "Aq", "NM", "UQ", "H0", "H1"),
                      what="isize")
gal6 <- readBamGappedAlignments(bamfile, param=param)
mcols(gal6)  # "tag" cols always after "what" cols

## ---------------------------------------------------------------------
## B. readBamGappedReads()
## ---------------------------------------------------------------------
greads1 <- readBamGappedReads(bamfile)
greads1
names(greads1)
qseq(greads1)
greads2 <- readBamGappedReads(bamfile, use.names=TRUE)
head(greads2)
head(names(greads2))

## ---------------------------------------------------------------------
## C. readBamGappedAlignmentPairs()
## ---------------------------------------------------------------------
galp1 <- readBamGappedAlignmentPairs(bamfile)
head(galp1)
names(galp1)
galp2 <- readBamGappedAlignmentPairs(bamfile, use.names=TRUE)
galp2
head(galp2)
head(names(galp2))
}

\keyword{manip}
