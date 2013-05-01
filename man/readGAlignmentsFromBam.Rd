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
\alias{groupAsPairs}
\alias{groupAsPairs,GAlignments-method}
\alias{groupAsPairs,GAlignmentsList-method}

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
readGAlignmentsFromBam(file, index=file, ..., use.names=FALSE, param=NULL)

readGappedReadsFromBam(file, index=file, use.names=FALSE, param=NULL)

readGAlignmentPairsFromBam(file, index=file, use.names=FALSE, param=NULL)

readGAlignmentsListFromBam(file, index=file, ..., use.names=FALSE,
                           param=ScanBamParam(), group.as.pairs=TRUE)

groupAsPairs(x, ...)
}

\arguments{
  \item{file, index}{The path to the BAM file to read, and to the index
    file of the BAM file to read, respectively. The latter is given
    \emph{without} the '.bai' extension. See \code{\link{scanBam}} for more
    information.}

  \item{x}{A \code{GAlignments} or \code{GAlignmentsList} object.
    The \code{GAlignments} must have names and metadata columns of
    \sQuote{flag}, \sQuote{mrnm} and \sQuote{mpos}.} 

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
    \code{readGAlignmentsFromBam} and \code{readGappedReadsFromBam}
    (i.e. only records corresponding to mapped reads are loaded),
    and \code{scanBamFlag(isUnmappedQuery=FALSE, isPaired=TRUE,
    hasUnmappedMate=FALSE)} for \code{readGAlignmentPairsFromBam}
    (i.e. only records corresponding to paired-end reads with both ends
    mapped are loaded).} No flag is set for \code{readGAlignmentsListFromBam}.

  \item{group.as.pairs}{A logical indicating if records should be 
    paired using the algorithm in \code{readGAlignmentPairsFromBam}.
    When \code{TRUE}, a metadata column is added to the result to
    indicate which records were paired using the algorithm. All other 
    (i.e., non-paired) records are grouped by QNAME.}
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
  and \code{?\link{makeGAlignmentPairs}} for the details of the
  pairing procedure.

  \code{readGAlignmentsListFromBam} returns a \code{GAlignmentsList} object. 
  When \code{group.as.pairs} is \code{TRUE} records are paired using the 
  algorithm in \code{readGAlignmentPairsFromBam}. A metadata column 
  is added to the result to indicate which records were paired using the 
  algorithm. Reads not paired by the algorithm are grouped by id 
  (QNAME in SAM/BAM file). When \code{group.as.pairs} is \code{FALSE}, no
  pairing is performed and all records are grouped by QNAME.
 
  See \code{?\link[GenomicRanges]{GAlignmentsList-class}} for a 
  description of \link[GenomicRanges]{GAlignmentsList} objects.

  \code{groupAsPairs} returns a \code{GAlignmentsList} object. 
  This function pairs reads in a \code{GAlignments} or 
  \code{GAlignmentsList} object using the algorithm in 
  \code{readGAlignmentPairsFromBam}. This function is called from
  \code{readGAlignmentsListFromBam} when \code{group.as.pairs} is
  TRUE.
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
  \code{readGAlignmentsListFromBam}. Single or paired-end data can be
  read into this structure. The list elements are groups of records
  by read id.
}

\note{
  BAM records corresponding to unmapped reads are always ignored.

  Starting with Rsamtools 1.7.1 (BioC 2.10), PCR or optical duplicates
  are loaded by default (use \code{scanBamFlag(isDuplicate=FALSE)} to
  drop them).
}

\author{H. Pages}

\seealso{
  \link[GenomicRanges]{GAlignments-class},
  \link[GenomicRanges]{GAlignmentsList-class},
  \link[ShortRead]{GappedReads-class},
  \link[GenomicRanges]{GAlignmentPairs-class},
  \code{\link{makeGAlignmentPairs}},
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

## Using the 'param' arg to drop PCR or optical duplicates and load
## additional BAM fields:
param <- ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE),
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
galp2 <- readGAlignmentPairsFromBam(bamfile, use.names=TRUE)
galp2
head(galp2)
head(names(galp2))

## ---------------------------------------------------------------------
## D. readGAlignmentsListFromBam()
## ---------------------------------------------------------------------

## (i) Records grouped by id (QUNAME):
## When 'yieldSize' is used, the Bam file must be sorted by qname.
## When creating the BamFile 'obeyQname' should be TRUE.
library(pasillaBamSubset)
fl <- untreated3_chr4() 
flsort <- sortBam(fl, tempfile(), byQname=TRUE)
bf <- BamFile(flsort, index=character(0), yieldSize=100, obeyQname=TRUE)
galist1 <- readGAlignmentsListFromBam(bf, group.as.pairs=FALSE)

## List elements are grouped by id and contain any number
## of pairs or fragments.
galist1
table(elementLengths(galist1))

## (ii) Records grouped with algorithm in readGAlignmentPairs():
fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
galist2 <- readGAlignmentsListFromBam(fl, group.as.pairs=TRUE)
## 'paired=TRUE' for records paired with the algorithm.
head(galist2, 2)
## Singletons or fragments not paired with the algorithm are grouped 
## by id.
tail(galist2, 2)

## groupAsPaired() applies the readGAlignmentsPairs() algorithm 
## to reads in a GAlignments or GAlignmentsList. These
## objects must have names and the appropriate metadata columns.
param <- ScanBamParam(what=c("flag", "mrnm", "mpos"))
gal <- readGAlignmentsFromBam(fl, use.names=TRUE, param=param) 
gal_paired <- groupAsPairs(gal)
}

\keyword{manip}
