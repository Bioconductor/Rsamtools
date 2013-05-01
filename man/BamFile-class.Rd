\name{BamFile}
\Rdversion{1.1}
\docType{class}
\alias{BamFile-class}
\alias{BamFileList-class}
% con/destructors
\alias{BamFile}
\alias{BamFileList}
\alias{open.BamFile}
\alias{close.BamFile}
% accessors
\alias{isOpen,BamFile-method}
% methods
\alias{scanBamHeader,BamFile-method}
\alias{seqinfo,BamFile-method}
\alias{obeyQname}
\alias{obeyQname<-}
\alias{obeyQname,BamFile-method}
\alias{obeyQname<-,BamFile-method}
\alias{obeyQname,BamFileList-method}
\alias{obeyQname<-,BamFileList-method}
\alias{scanBam,BamFile-method}
\alias{countBam,BamFile-method}
\alias{countBam,BamFileList-method}
\alias{filterBam,BamFile-method}
\alias{indexBam,BamFile-method}
\alias{sortBam,BamFile-method}
\alias{readGAlignmentsFromBam,BamFile-method}
\alias{readGappedReadsFromBam,BamFile-method}
\alias{readGAlignmentPairsFromBam,BamFile-method}
\alias{readGAlignmentsListFromBam,BamFile-method}

\alias{mergeBam,BamFileList-method}
% counting 
\alias{summarizeOverlaps,GRanges,BamFileList-method}
\alias{summarizeOverlaps,GRangesList,BamFileList-method}
\alias{findSpliceOverlaps,character-method}
\alias{findSpliceOverlaps,BamFile-method}
\alias{findSpliceOverlaps,character,ANY-method}
\alias{findSpliceOverlaps,BamFile,ANY-method}
\alias{quickCountBam,BamFile-method}
\alias{coverage,BamFile-method}

\alias{show,BamFile-method}
\alias{show,BamFileList-method}

\title{Maintain and use BAM files}

\description{

  Use \code{BamFile()} to create a reference to a BAM file (and
  optionally its index). The reference remains open across calls to
  methods, avoiding costly index re-loading.

  \code{BamFileList()} provides a convenient way of managing a list of
  \code{BamFile} instances.

}

\usage{

## Constructors

BamFile(file, index=file, ..., yieldSize=NA_integer_, obeyQname=FALSE)
BamFileList(..., yieldSize=NA_integer_, obeyQname=FALSE)

## Opening / closing

\S3method{open}{BamFile}(con, ...)
\S3method{close}{BamFile}(con, ...)

## accessors; also path(), index(), yieldSize()

\S4method{isOpen}{BamFile}(con, rw="")
\S4method{obeyQname}{BamFile}(object, ...)
obeyQname(object, ...) <- value

## actions

\S4method{scanBamHeader}{BamFile}(files, ...)
\S4method{seqinfo}{BamFile}(x)
\S4method{scanBam}{BamFile}(file, index=file, ..., param=ScanBamParam(what=scanBamWhat()))
\S4method{countBam}{BamFile}(file, index=file, ..., param=ScanBamParam())
\S4method{countBam}{BamFileList}(file, index=file, ..., param=ScanBamParam())
\S4method{filterBam}{BamFile}(file, destination, index=file, ...,
    indexDestination=TRUE, param=ScanBamParam(what=scanBamWhat()))
\S4method{indexBam}{BamFile}(files, ...)
\S4method{sortBam}{BamFile}(file, destination, ..., byQname=FALSE, maxMemory=512)
\S4method{mergeBam}{BamFileList}(files, destination, ...)
\S4method{readGAlignmentsFromBam}{BamFile}(file, index=file, ..., use.names=FALSE, param=NULL)
\S4method{readGappedReadsFromBam}{BamFile}(file, index=file, use.names=FALSE, param=NULL)
\S4method{readGAlignmentPairsFromBam}{BamFile}(file, index=file, use.names=FALSE, param=NULL)
\S4method{readGAlignmentsListFromBam}{BamFile}(file, index=file, ...,
use.names=FALSE, param=ScanBamParam(), group.as.pairs=TRUE)

## counting

\S4method{summarizeOverlaps}{GRanges,BamFileList}(features, reads, mode, ignore.strand=FALSE, ...,
    singleEnd=TRUE, param=ScanBamParam()) 
\S4method{summarizeOverlaps}{GRangesList,BamFileList}(features, reads, mode, ignore.strand=FALSE, ...,
    singleEnd=TRUE, param=ScanBamParam()) 

\S4method{findSpliceOverlaps}{character,ANY}(query, subject, ignore.strand=FALSE, ...,
    param=ScanBamParam(), pairedEnd=FALSE)
\S4method{findSpliceOverlaps}{BamFile,ANY}(query, subject, ignore.strand=FALSE, ...,
    param=ScanBamParam(), pairedEnd=FALSE)

\S4method{coverage}{BamFile}(x, shift=0L, width=NULL, weight=1L, ..., param = ScanBamParam())

\S4method{quickCountBam}{BamFile}(file, ..., param=ScanBamParam(), mainGroupsOnly=FALSE)
}

\arguments{

  \item{...}{Additional arguments. For \code{BamFileList}, this can
    either be a single character vector of paths to BAM files, or
    several instances of \code{BamFile} objects. When a character vector
    of paths, a second named argument \sQuote{index} can be a
    \code{character()} vector of length equal to the first argument
    specifying the paths to the index files, or \code{character()} to
    indicate that no index file is available. See
    \code{\link{BamFile}}. For \code{coverage}, the arguments are passed
    to the
    \code{\link[GenomicRanges:coverage,GAlignments-method]{coverage}}
    method for \code{GAlignments} objects.}

  \item{con}{An instance of \code{BamFile}.}

  \item{x, object, file, files}{A character vector of BAM file paths
    (for \code{BamFile}) or a \code{BamFile} instance (for other
    methods).}

  \item{index}{character(1); the BAM index file path (for
    \code{BamFile}); ignored for all other methods on this page.}

  \item{yieldSize}{Number of records to yield each time the file is read
    from using \code{scanBam}. Only valid when
    \code{length(bamWhich(param)) == 0}. \code{yieldSize} does not alter
    existing yield sizes, include \code{NA}, when creating a
    \code{BamFileList} from \code{BamFile} instances.}

  \item{destination}{character(1) file path to write filtered reads to.}

  \item{indexDestination}{logical(1) indicating whether the destination
    file should also be indexed.}

  \item{byQname, maxMemory}{See \code{\link{sortBam}}.}

  \item{obeyQname, value}{A logical(1) indicating whether the file is
    sorted by \code{qname}.}

  \item{param}{An optional \code{\linkS4class{ScanBamParam}} instance to
     further influence scanning, counting, or filtering.}

  \item{use.names}{Construct the names of the returned object
     from the query template names (QNAME field)? If not (the default),
     then the returned object has no names.}

  \item{rw}{Mode of file; ignored.}

  \item{reads}{
    A \code{\linkS4class{BamFileList}} that represents the data to be
    counted by \code{summarizeOverlaps}.}

  \item{features}{
    A \link{GRanges} or a \link{GRangesList} object of genomic regions of 
    interest. When a \link{GRanges} is supplied, each row is considered a 
    feature. When a \link{GRangesList} is supplied, each higher list-level is
    considered a feature. This distinction is important when defining an overlap
    between a read and a feature. See ?\code{summarizeOverlaps} for details.}

  \item{mode}{
    A function that defines the method to be used when a read overlaps
    more than one feature. Pre-defined options are "Union",
    "IntersectionStrict", or "IntersectionNotEmpty" and are designed
    after the counting modes available in the HTSeq package by Simon
    Anders (see references).

    \itemize{
      \item "Union" : (Default) Reads that overlap any portion of exactly one 
            feature are counted. Reads that overlap multiple features are 
            discarded. 
      \item "IntersectionStrict" : A read must fall completely "within" the
            feature to be counted. If a read overlaps multiple features but
            falls "within" only one, the read is counted for that feature.
            If the read is "within" multiple features, the read is discarded. 
      \item "IntersectionNotEmpty" : A read must fall in a unique disjoint
            region of a feature to be counted. When a read overlaps multiple
            features, the features are partitioned into disjoint intervals. 
            Regions that are shared between the features are discarded leaving
            only the unique disjoint regions. If the read overlaps one of 
            these remaining regions, it is assigned to the feature the
            unique disjoint region came from.
   }}

  \item{ignore.strand}{
    A logical value indicating if strand should be considered when matching.}

  \item{singleEnd}{
    A logical value indicating if reads are single or paired-end.}

  \item{pairedEnd}{
    A logical value indicating if reads are single or paired-end.}

  \item{query}{
    \code{character} name of a Bam file, a \link[Rsamtools]{BamFile},
    \linkS4class{GAlignments}, \linkS4class{GAlignmentPairs}
    or a \linkS4class{GRangesList} object containing the reads.

    Paired-end reads can be supplied in a Bam file or 
    \linkS4class{GAlignmentPairs} object. Single-end may
    be in a Bam file, \linkS4class{GAlignments} or \linkS4class{GRanges}
    object.}

  \item{subject}{
    A \link[GenomicFeatures]{TranscriptDb}, or \link{GRangesList} containing
    the annotations.}

  \item{shift, width, weight}{See \code{\link[IRanges:coverage]{coverage}}.}

  \item{mainGroupsOnly}{See \code{\link{quickCountBam}}.}

  \item{group.as.pairs}{A logical indicating if the records should be paired
    using the algorithm in \code{readGAlignmentPairs}. When TRUE, a metadata 
    column named \code{paired} is returned with the result indicating
    which records were paired using the algorithm. All other reads are 
    grouped by read id (QNAME in SAM/BAM). Applies to 
    \code{readGAlignmentsListFromBam} only.}
}

\section{Objects from the Class}{

  Objects are created by calls of the form \code{BamFile()}.

}

\section{Fields}{

  The \code{BamFile} class inherits fields from the
  \code{\linkS4class{RsamtoolsFile}} class and has field:

  \describe{

    \item{obeyQname}{A logical(0) indicating if the file was sorted by qname.}

  }
}

\section{Functions and methods}{

  \code{BamFileList} inherits methods from
  \code{\link{RsamtoolsFileList}} and \code{\link{SimpleList}}.

  Opening / closing:
  \describe{

    \item{open.BamFile}{Opens the (local or remote) \code{path} and
      \code{index} (if \code{bamIndex} is not \code{character(0)}),
      files.  Returns a \code{BamFile} instance.}

    \item{close.BamFile}{Closes the \code{BamFile} \code{con}; returning
      (invisibly) the updated \code{BamFile}. The instance may be
      re-opened with \code{open.BamFile}.}

  }

  Accessors: 
  \describe{

    \item{path}{Returns a character(1) vector of BAM path names.}

    \item{index}{Returns a character(1) vector of BAM index path
      names.}

    \item{yieldSize, yieldSize<-}{Return or set an integer(1) vector
      indicating yield size.}

    \item{obeyQname, obeyQname<-}{Return or set a logical(0)
      indicating if the file was sorted by qname.}
  }

  Methods:
  \describe{

    \item{scanBamHeader}{Visit the path in \code{path(file)}, returning
      the information contained in the file header; see
      \code{\link{scanBamHeader}}.}

    \item{seqinfo}{Visit the path in \code{path(file)}, returning
      a \code{\linkS4class{Seqinfo}} instance containing information on
      the lengths of each sequence.}

    \item{scanBam}{Visit the path in \code{path(file)}, returning the
      result of \code{\link{scanBam}} applied to the specified path.}

    \item{countBam}{Visit the path(s) in \code{path(file)}, returning
      the result of \code{\link{countBam}} applied to the specified
      path.}

    \item{filterBam}{Visit the path in \code{path(file)}, returning
      the result of \code{\link{filterBam}} applied to the specified
      path.}

    \item{indexBam}{Visit the path in \code{path(file)}, returning
      the result of \code{\link{indexBam}} applied to the specified
      path.}

    \item{sortBam}{Visit the path in \code{path(file)}, returning the
      result of \code{\link{sortBam}} applied to the specified path.}

    \item{mergeBam}{Merge several BAM files into a single BAM file. See
      \code{\link{mergeBam}} for details; additional arguments supported
      by \code{mergeBam,character-method} are also available for
      \code{BamFileList}.}

    \item{readGAlignmentsFromBam, readGappedReadsFromBam,
          readGAlignmentPairsFromBam}{
      Visit the path in \code{path(file)}, returning the result of
      \code{readGAlignmentsFromBam}, \code{readGappedReadsFromBam},
      or \code{readGAlignmentPairsFromBam} applied to the specified path.
      See \code{\link{readGAlignmentsFromBam}}.}

    \item{readGAlignmentsListFromBam}{
      Visit the Bam file in \code{path(file)}. When \code{yieldSize} is
      set on the BamFile, the file must be sorted by qname, 
      see ?\code{sortBam}. For a full description of the method see
      \code{\link{readGAlignmentsListFromBam}}.} 

    \item{show}{Compactly display the object.}

  }

}
\author{Martin Morgan and Marc Carlson}

\seealso{
  The \code{GenomicRanges} package is where the \code{summarizeOverlaps}
  method originates.
}

\examples{

fl <- system.file("extdata", "ex1.bam", package="Rsamtools",
                  mustWork=TRUE)
length(scanBam(fl)[[1]][[1]])  # all records

bf <- open(BamFile(fl))        # implicit index
bf
identical(scanBam(bf), scanBam(fl))
close(bf)

## chunks of size 1000
bf <- open(BamFile(fl, yieldSize=1000)) 
while (nrec <- length(scanBam(bf)[[1]][[1]]))
    cat("records:", nrec, "\n")
close(bf)

rng <- GRanges(c("seq1", "seq2"), IRanges(1, c(1575, 1584)))

## repeatedly visit 'bf'
bf <- open(BamFile(fl))
sapply(seq_len(length(rng)), function(i, bamFile, rng) {
    param <- ScanBamParam(which=rng[i], what="seq")
    bam <- scanBam(bamFile, param=param)[[1]]
    alphabetFrequency(bam[["seq"]], baseOnly=TRUE, collapse=TRUE)
}, bf, rng)
close(bf)


##------------------------------------------------------------------------
## summarizeOverlaps with BamFileList
##

library(pasillaBamSubset)
library("TxDb.Dmelanogaster.UCSC.dm3.ensGene")
exbygene <- exonsBy(TxDb.Dmelanogaster.UCSC.dm3.ensGene, "gene")

## single-end:
## When 'yieldSize' is specified the file is processed by chunks. 
## Otherwise the complete file is read into memory.
fl <- untreated1_chr4() 
bfl <- BamFileList(fl, yieldSize=50000)
se1 <- summarizeOverlaps(exbygene, bfl, singleEnd=TRUE)
counts1 <- assays(se1)$counts

## paired-end sorted by qname:
## Set 'singleEnd' to 'FALSE'. A BAM file sorted by qname
## can be read in chunks with 'yieldSize'. 
fl <- untreated3_chr4() 
sortfl <- sortBam(fl, tempfile(), byQname=TRUE)
bf2 <- BamFileList(sortfl, index=character(0), 
                   yieldSize=50000, obeyQname=TRUE)
se2 <- summarizeOverlaps(exbygene, bf2, singleEnd=FALSE)
counts2 <- assays(se2)$counts

## paired-end not sorted:
## If the file is not sorted by qname, all records are read 
## into memory for sorting and to determine proper pairs.
## Any 'yieldSize' set on the BamFile will be ignored.
fl <- untreated3_chr4() 
bf3 <- BamFileList(fl)
se3 <- summarizeOverlaps(exbygene, bf3, singleEnd=FALSE)
counts3 <- assays(se3)$counts

identical(as.vector(counts2), as.vector(counts3))

##------------------------------------------------------------------------
## findSpliceOverlaps 
##

## See ?'findSpliceOverlaps' for examples 
}

\keyword{classes}
