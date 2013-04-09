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
readGAlignmentsFromBam(file, index=file, ..., use.names=FALSE, param=NULL)

readGappedReadsFromBam(file, index=file, use.names=FALSE, param=NULL)

readGAlignmentPairsFromBam(file, index=file, use.names=FALSE, param=NULL)

readGAlignmentsListFromBam(file, index=file, ..., use.names=FALSE,
                           param=ScanBamParam(), asProperPairs=TRUE)
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
    \code{readGAlignmentsFromBam} and \code{readGappedReadsFromBam}
    (i.e. only records corresponding to mapped reads are loaded),
    and \code{scanBamFlag(isUnmappedQuery=FALSE, isPaired=TRUE,
    hasUnmappedMate=FALSE)} for \code{readGAlignmentPairsFromBam}
    (i.e. only records corresponding to paired-end reads with both ends
    mapped are loaded).} No flag is set for \code{readGAlignmentsListFromBam}
    unless \code{asProperPairs=TRUE}. In this case, the records read in
    are the same as for \code{readGAlignmentPairsFromBam}.

  \item{asProperPairs}{A logical indicating if the records should be filtered
    such that only proper pairs are returned. Applies to 
    \code{readGAlignmentsListFromBam} only. If \code{TRUE}, the records 
    are the same as from a call to \code{readGAlignmentPairsFromBam} 
    except a \code{GAlignmentsList} is returned instead of a 
    \code{GAlignmentPairs} object.}
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

  The return value from \code{readGAlignmentsListFromBam} is a 
  \code{GAlignmentsList} where each list element contains all records 
  of the same id (QNAME in SAM/BAM file). When \code{asProperPairs} is 
  \code{TRUE} each list element has exactly 2 records; these are the 
  same data as that returned from \code{readGAlignmentPairsFromBam}, only 
  the return class is different. When \code{asProperPairs} is \code{FALSE}, 
  no QC is performed resulting in 1 or more records per element. List 
  elements containing singletons, unpaired reads or single fragments have 
  a length of 1 while paired-end reads or those with multiple fragments 
  have a length of 2 or greater.
  (NOTE: asProperPairs=TRUE not yet implemented)

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
## D. readGAlignmentPairsFromBam()
## ---------------------------------------------------------------------
## As sample data we use the paired-end file from pasillaBamSubset.
## This method requires that Bam file to be sorted by qname. Setting
## the yieldSize is optional.
library(pasillaBamSubset)
bfsort <- sortBam(untreated3_chr4(), tempfile(), byQname=TRUE)
bf <- BamFile(bfsort, index=character(0), yieldSize=100, obeyQname=TRUE)
galist <- readGAlignmentsListFromBam(bf, asProperPairs=FALSE)

## List elements are grouped by id and can hold any number
## of pairs or fragments.
galist
table(elementLengths(galist))

## Single reads appearing in 'galist' are filtered out by 
## the readGAlignmentPairsFromBam QC process because they
## are not proper pairs. 
galp <- readGAlignmentPairsFromBam(bf)
findOverlaps(galist[elementLengths(galist) == 1], 
             unlist(galp), type="within")

## When 'asProperPairs' is 'TRUE', readGAlignmentPairsFromBam 
## and readGAlignmentsList return the same records but in different
## classes.
}

\keyword{manip}
