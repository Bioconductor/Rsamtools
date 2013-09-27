\name{readGAlignmentsFromBam}
\Rdversion{1.1}

\alias{readGAlignmentsFromBam}
\alias{readGAlignmentsFromBam,character-method}
\alias{readGappedReadsFromBam}
\alias{readGappedReadsFromBam,character-method}
\alias{readGAlignmentPairsFromBam}
\alias{readGAlignmentPairsFromBam,character-method}
\alias{readGAlignmentsListFromBam}
\alias{readGAlignmentsListFromBam,character-method}

% Old stuff:
\alias{readBamGappedAlignments}
\alias{readBamGappedReads}
\alias{readBamGappedAlignmentPairs}
\alias{readBamGAlignmentsList}

\title{Reading a GAlignments, GappedReads, GAlignmentPairs, or GAlignmentsList
       object from a BAM file}

\description{
Read a \link[GenomicRanges]{GAlignments}, \link[ShortRead]{GappedReads},
\link[GenomicRanges]{GAlignmentPairs}, or \link[GenomicRanges]{GAlignmentsList}
object from a BAM file.
}

\usage{
readGAlignmentsFromBam(file, index=file, ..., use.names=FALSE, param=NULL,
                       with.which_label=FALSE)

readGappedReadsFromBam(file, index=file, use.names=FALSE, param=NULL,
                       with.which_label=FALSE)

readGAlignmentPairsFromBam(file, index=file, use.names=FALSE, param=NULL,
                           with.which_label=FALSE)

readGAlignmentsListFromBam(file, index=file, ..., use.names=FALSE,
                           param=ScanBamParam())
}

\arguments{
  \item{file, index}{The path to the BAM file to read, and to the index
    file of the BAM file to read, respectively. The latter is given
    \emph{without} the '.bai' extension. See \code{\link{scanBam}} for more
    information.}

  \item{\dots}{Arguments passed to other methods.}

  \item{use.names}{Use the query template names (QNAME field) as the
    names of the returned object? If not (the default), then the returned
    object has no names.}

  \item{param}{\code{NULL} or an instance of \code{\linkS4class{ScanBamParam}}.
    Like for \code{\link{scanBam}}, this influences what fields and which
    records are imported. However, note that the fields specified
    thru this \code{\linkS4class{ScanBamParam}} object will be loaded
    \emph{in addition} to any field required for generating the returned
    object (\link[GenomicRanges]{GAlignments},
    \link[ShortRead]{GappedReads}, or
    \link[GenomicRanges]{GAlignmentPairs} object), but only the fields
    requested by the user will actually be kept as metadata columns of the
    object.

    By default (i.e. \code{param=NULL} or \code{param=ScanBamParam()}), no 
    additional field is loaded. The flag used is 
    \code{scanBamFlag(isUnmappedQuery=FALSE)} for
    \code{readGAlignmentsFromBam}, \code{readGappedReadsFromBam} and
    \code{readGAlignmentsListFromBam}
    (i.e. only records corresponding to mapped reads are loaded),
    and \code{scanBamFlag(isUnmappedQuery=FALSE, isPaired=TRUE,
    hasUnmappedMate=FALSE)} for \code{readGAlignmentPairsFromBam}
    (i.e. only records corresponding to paired-end reads with both ends
    mapped are loaded).}

  \item{with.which_label}{\code{TRUE} or \code{FALSE}. Only supported
    by \code{readGAlignmentsFromBam}, \code{readGappedReadsFromBam},
    and \code{readGAlignmentPairsFromBam} at the moment. If \code{TRUE}
    and if \code{param} has a \code{which} component, a \code{"which_label"}
    metadata column is added to the returned \link[GenomicRanges]{GAlignments}
    or \link[ShortRead]{GappedReads} object, or to the \code{\link{first}}
    and \code{\link{last}} components of the returned
    \link[GenomicRanges]{GAlignmentPairs} object.
    The purpose of this metadata column is to unambiguously identify
    the range in \code{which} where each element in the returned object
    originates from. The labels used to identify the ranges are normally
    of the form \code{"seq1:12250-246500"}, that is, they're the same as
    the names found on the outer list that \code{\link{scanBam}} would
    return if called with the same \code{param} argument. If some ranges
    are duplicated, then the labels are made unique by appending a unique
    suffix to all of them. The \code{"which_label"} metadata column is
    represented as a factor-\link[IRanges]{Rle}.}
}

\details{
  See \code{?\link[GenomicRanges]{GAlignments-class}} for a
  description of \link[GenomicRanges]{GAlignments} objects.

  See \code{?\link[ShortRead]{GappedReads-class}} for a
  description of \link[ShortRead]{GappedReads} objects.

  \code{readGAlignmentPairsFromBam} proceeds in 2 steps:
  \enumerate{
    \item Load the BAM file into a \link[GenomicRanges]{GAlignments}
          object with \code{readGAlignmentsFromBam};
    \item Turn this \link[GenomicRanges]{GAlignments} object into a
          \link[GenomicRanges]{GAlignmentPairs} object by pairing
          its elements.
  }
  See \code{?\link[GenomicRanges]{GAlignmentPairs-class}} for a
  description of \link[GenomicRanges]{GAlignmentPairs} objects,
  and \code{?\link{findMateAlignment}} for a description of the
  pairing algorithm (including timing and memory requirement).

  \code{readGAlignmentsListFromBam} pairs records into \sQuote{mates}
  acording to the criteria below. A \code{GAlignmentsList} is 
  returned with a metadata column, \sQuote{mates}, indicating mate 
  status. The mates are returned first followed by non-mates. When
  the \sQuote{file} argument is a BamFile, \sQuote{asMates=TRUE}
  must be set. See the \sQuote{asMates} section of 
  \code{?\link{BamFile}} for details. 

    Mate criteria:
    \itemize{
      \item {Bit 0x1 (multiple segments) is 1}.
      \item {Bit 0x4 (segment unmapped) is 0}.
      \item {Bit 0x8 (next segment unmapped) is 0}.
      \item {Bit 0x40 and 0x80 (first/last segment): 
             Segments are a pair of first/last OR
             neither segment is marked first/last}.
      \item {Bit 0x100 (secondary alignment): 
             Both segments are secondary OR both not secondary}
      \item {\sQuote{qname} match}.
      \item {\sQuote{tid} match}.
      \item {segment1 \sQuote{mpos} matches segment2 \sQuote{pos} AND
             segment2 \sQuote{mpos} matches segment1 \sQuote{pos}}
    } 
  Records that do not pass these criteria are returned with
  mate status FALSE. Flags, tags and ranges may be specified in the
  \code{ScanBamParam} for fine tuning of results.

  See \code{?\link[GenomicRanges]{GAlignmentsList-class}} for a 
  description of \link[GenomicRanges]{GAlignmentsList} objects.
}

\value{
  A \link[GenomicRanges]{GAlignments} object for
  \code{readGAlignmentsFromBam}.

  A \link[ShortRead]{GappedReads} object for \code{readGappedReadsFromBam}.

  A \link[GenomicRanges]{GAlignmentPairs} object for
  \code{readGAlignmentPairsFromBam}.
  Note that a BAM (or SAM) file can in theory contain a mix of single-end
  and paired-end reads, but in practise it seems that single-end and
  paired-end are not mixed. In other words, the value of flag bit 0x1
  (\code{isPaired}) is the same for all the records in a file.
  So if \code{readGAlignmentPairsFromBam} returns a
  \link[GenomicRanges]{GAlignmentPairs} object of length zero,
  this almost certainly means that the BAM (or SAM) file contains
  alignments for single-end reads (although it could also mean that the
  user-supplied \code{\linkS4class{ScanBamParam}} is filtering out everything,
  or that the file is empty, or that all the records in the file correspond
  to unmapped reads).

  A \link[GenomicRanges]{GAlignmentsList} object for
  \code{readGAlignmentsListFromBam}. 
}

\note{
  BAM records corresponding to unmapped reads are always ignored.

  Starting with Rsamtools 1.7.1 (BioC 2.10), PCR or optical duplicates
  are loaded by default (use \code{scanBamFlag(isDuplicate=FALSE)} to
  drop them).
}

\author{H. Pages <hpages@fhcrc.org> and
        Valerie Obenchain <vobencha@fhcrc.org>}

\seealso{
  \link[GenomicRanges]{GAlignments-class},
  \link[GenomicRanges]{GAlignmentsList-class},
  \link[ShortRead]{GappedReads-class},
  \link[GenomicRanges]{GAlignmentPairs-class},
  \code{\link{findMateAlignment}},
  \code{\link{scanBam}},
  \code{\link{ScanBamParam}}
}

\examples{
## ---------------------------------------------------------------------
## A. readGAlignmentsFromBam()
## ---------------------------------------------------------------------

## Simple use:
bamfile <- system.file("extdata", "ex1.bam", package="Rsamtools",
                       mustWork=TRUE)
gal1 <- readGAlignmentsFromBam(bamfile)
gal1
names(gal1)

## Using the 'use.names' arg:
gal2 <- readGAlignmentsFromBam(bamfile, use.names=TRUE)
gal2
head(names(gal2))

## Using the 'param' arg to drop PCR or optical duplicates as well as
## secondary alignments, and to load additional BAM fields:
param <- ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE,
                                       isNotPrimaryRead=FALSE),
                      what=c("qual", "flag"))
gal3 <- readGAlignmentsFromBam(bamfile, param=param)
gal3
mcols(gal3)

## Using the 'param' arg to load reads from particular regions.
## Note that if we weren't providing a 'what' argument here, all the
## BAM fields would be loaded:
which <- RangesList(seq1=IRanges(1000, 2000),
                    seq2=IRanges(c(100, 1000), c(1000, 2000)))
param <- ScanBamParam(which=which)
gal4 <- readGAlignmentsFromBam(bamfile, param=param)
gal4

## Note that a given record is loaded one time for each region it
## belongs to (this is a scanBam() feature, readGAlignmentsFromBam()
## is based on scanBam()):
which <- IRangesList(seq2=IRanges(c(1563, 1567), width=1))
param <- ScanBamParam(which=which)
gal5 <- readGAlignmentsFromBam(bamfile, param=param)
gal5

## Use 'with.which_label=TRUE' to identify the range in 'which'
## where each element in 'gal5' originates from.
gal5 <- readGAlignmentsFromBam(bamfile, param=param,
                               with.which_label=TRUE)
gal5

## Using the 'param' arg to load tags. Except for MF and Aq, the tags
## specified below are predefined tags (see the SAM Spec for the list
## of predefined tags and their meaning).
param <- ScanBamParam(tag=c("MF", "Aq", "NM", "UQ", "H0", "H1"),
                      what="isize")
gal6 <- readGAlignmentsFromBam(bamfile, param=param)
mcols(gal6)  # "tag" cols always after "what" cols

## ---------------------------------------------------------------------
## B. readGappedReadsFromBam()
## ---------------------------------------------------------------------
greads1 <- readGappedReadsFromBam(bamfile)
greads1
names(greads1)
qseq(greads1)
greads2 <- readGappedReadsFromBam(bamfile, use.names=TRUE)
head(greads2)
head(names(greads2))

## ---------------------------------------------------------------------
## C. readGAlignmentPairsFromBam()
## ---------------------------------------------------------------------
galp1 <- readGAlignmentPairsFromBam(bamfile)
head(galp1)
names(galp1)
## Using the 'param' arg to drop PCR or optical duplicates as well as
## secondary alignments (dropping secondary alignments can help make the
## pairing algorithm run significantly faster, see ?findMateAlignment):
param <- ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE,
                                       isNotPrimaryRead=FALSE))
galp2 <- readGAlignmentPairsFromBam(bamfile, use.names=TRUE, param=param)
galp2
head(galp2)
head(names(galp2))

## ---------------------------------------------------------------------
## D. readGAlignmentsListFromBam()
## ---------------------------------------------------------------------

library(pasillaBamSubset)

## 'file' as character.
fl <- untreated3_chr4() 
galist1 <- readGAlignmentsListFromBam(fl)
galist1[1:3]
length(galist1)
table(elementLengths(galist1))

## When 'file' is a BamFile, 'asMates' must be TRUE. If FALSE,
## the data are treated as single-end and each list element of the
## GAlignmentsList will be of length 1. For single-end data 
## use readGAlignements().
bf <- BamFile(fl, yieldSize=3, asMates=TRUE)
readGAlignmentsList(bf)

## Use a 'param' to fine tune the results.
param <- ScanBamParam(flag=scanBamFlag(isProperPair=TRUE))
galist2 <- readGAlignmentsListFromBam(fl, param=param)
length(galist2)
}

\keyword{manip}
