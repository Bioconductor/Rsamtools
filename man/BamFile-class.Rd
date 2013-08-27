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
\alias{isIncomplete,BamFile-method}
% methods
\alias{scanBamHeader,BamFile-method}
\alias{seqinfo,BamFile-method}
\alias{obeyQname}
\alias{obeyQname<-}
\alias{obeyQname,BamFile-method}
\alias{obeyQname<-,BamFile-method}
\alias{obeyQname,BamFileList-method}
\alias{obeyQname<-,BamFileList-method}
\alias{asMates}
\alias{asMates<-}
\alias{asMates,BamFile-method}
\alias{asMates<-,BamFile-method}
\alias{asMates,BamFileList-method}
\alias{asMates<-,BamFileList-method}
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
\alias{summarizeOverlaps,GRanges,character-method}
\alias{summarizeOverlaps,GRangesList,character-method}
\alias{summarizeOverlaps,GRanges,BamFileList-method}
\alias{summarizeOverlaps,GRangesList,BamFileList-method}
\alias{summarizeOverlaps,GRanges,BamFile-method}
\alias{summarizeOverlaps,GRangesList,BamFile-method}
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

BamFile(file, index=file, ..., yieldSize=NA_integer_, obeyQname=FALSE,
        asMates=FALSE)
BamFileList(..., yieldSize=NA_integer_, obeyQname=FALSE, asMates=FALSE)

## Opening / closing

\S3method{open}{BamFile}(con, ...)
\S3method{close}{BamFile}(con, ...)

## accessors; also path(), index(), yieldSize()

\S4method{isOpen}{BamFile}(con, rw="")
\S4method{isIncomplete}{BamFile}(con)
\S4method{obeyQname}{BamFile}(object, ...)
obeyQname(object, ...) <- value
\S4method{asMates}{BamFile}(object, ...)
asMates(object, ...) <- value

## actions

\S4method{scanBamHeader}{BamFile}(files, ...)
\S4method{seqinfo}{BamFile}(x)
\S4method{scanBam}{BamFile}(file, index=file, ..., param=ScanBamParam(what=scanBamWhat()))
\S4method{countBam}{BamFile}(file, index=file, ..., param=ScanBamParam())
\S4method{countBam}{BamFileList}(file, index=file, ..., param=ScanBamParam())
\S4method{filterBam}{BamFile}(file, destination, index=file, ...,
    filter=FilterRules(), indexDestination=TRUE,
    param=ScanBamParam(what=scanBamWhat()))
\S4method{indexBam}{BamFile}(files, ...)
\S4method{sortBam}{BamFile}(file, destination, ..., byQname=FALSE, maxMemory=512)
\S4method{mergeBam}{BamFileList}(files, destination, ...)
\S4method{readGAlignmentsFromBam}{BamFile}(file, index=file, ..., use.names=FALSE, param=NULL)
\S4method{readGappedReadsFromBam}{BamFile}(file, index=file, use.names=FALSE, param=NULL)
\S4method{readGAlignmentPairsFromBam}{BamFile}(file, index=file, use.names=FALSE, param=NULL)
\S4method{readGAlignmentsListFromBam}{BamFile}(file, index=file, ...,
use.names=FALSE, param=ScanBamParam(), group.as.pairs=TRUE)

## counting

\S4method{summarizeOverlaps}{GRanges,character}(features, reads, mode,
    ignore.strand=FALSE, ..., yieldSize=1000000L, inter.feature=TRUE,
    singleEnd=TRUE, fragments=TRUE, param=ScanBamParam())
\S4method{summarizeOverlaps}{GRangesList,character}(features, reads, mode,
    ignore.strand=FALSE, ..., yieldSize=1000000L, inter.feature=TRUE,
    singleEnd=TRUE, fragments=TRUE, param=ScanBamParam())
\S4method{summarizeOverlaps}{GRanges,BamFileList}(features, reads, mode, 
    ignore.strand=FALSE, ..., inter.feature=TRUE, singleEnd=TRUE, 
    fragments=TRUE, param=ScanBamParam()) 
\S4method{summarizeOverlaps}{GRangesList,BamFileList}(features, reads, mode, 
    ignore.strand=FALSE, ..., inter.feature=TRUE, singleEnd=TRUE, 
    fragments=TRUE, param=ScanBamParam()) 
\S4method{summarizeOverlaps}{GRanges,BamFile}(features, reads, mode, 
    ignore.strand=FALSE, ..., inter.feature=TRUE, singleEnd=TRUE, 
    fragments=TRUE, param=ScanBamParam()) 
\S4method{summarizeOverlaps}{GRangesList,BamFile}(features, reads, mode, 
    ignore.strand=FALSE, ..., inter.feature=TRUE, singleEnd=TRUE, 
    fragments=TRUE, param=ScanBamParam()) 

\S4method{findSpliceOverlaps}{character,ANY}(query, subject, ignore.strand=FALSE, ...,
    param=ScanBamParam(), singleEnd=TRUE)
\S4method{findSpliceOverlaps}{BamFile,ANY}(query, subject, ignore.strand=FALSE, ...,
    param=ScanBamParam(), singleEnd=TRUE)

\S4method{coverage}{BamFile}(x, shift=0L, width=NULL, weight=1L, ..., param = ScanBamParam())

\S4method{quickCountBam}{BamFile}(file, ..., param=ScanBamParam(), mainGroupsOnly=FALSE)
}

\arguments{

  \item{...}{Additional arguments.

    For \code{BamFileList}, this can either be a single character vector
    of paths to BAM files, or several instances of \code{BamFile}
    objects. When a character vector of paths, a second named argument
    \sQuote{index} can be a \code{character()} vector of length equal to
    the first argument specifying the paths to the index files, or
    \code{character()} to indicate that no index file is available. See
    \code{\link{BamFile}}.

    For \code{coverage}, the arguments are passed to the
    \code{\link[GenomicRanges:coverage,GAlignments-method]{coverage}}
    method for \code{GAlignments} objects.

    For \code{summarizeOverlaps}, providing
    \code{count.mapped.reads=TRUE} include additional passes through the
    BAM file to collect statistics like those from \code{countBam}.

  }

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

  \item{filter}{A \code{\link{FilterRules}} instance. Functions in the
    \code{FilterRules} instance should expect a single \code{DataFrame}
    argument representing all information specified by
    \code{param}. Each function must return a \code{logical} vector,
    usually of length equal to the number of rows of the
    \code{DataFrame}. Return values are used to include (when
    \code{TRUE}) corresponding records in the filtered BAM file.}

  \item{destination}{character(1) file path to write filtered reads to.}

  \item{indexDestination}{logical(1) indicating whether the destination
    file should also be indexed.}

  \item{byQname, maxMemory}{See \code{\link{sortBam}}.}

  \item{obeyQname, value}{A logical(1) indicating whether the file is
    sorted by \code{qname}.}

  \item{asMates}{A logical(1) indicating whether records should
    be returned as mated pairs. When \code{TRUE} scanBam attempts
    to pair, or \sQuote{mate}, the records and returns two additional 
    fields of \code{partition} and \code{mates}. The \code{partition}
    vector indicates the number of records in each group, similar to
    \sQuote{width} in a partition object. The \code{mates} vector
    is 1 for records that have been mated by the algorithm and 0 for
    non-mates. The non-mates are grouped by \sQuote{qname} and may
    be returned in groups of 1 or more (i.e., the \code{partition}
    may be greater than 1 but the \code{mates} value will be 0 for all).

    Records are considered \sQuote{mates} (i.e., mate status 1) if
    they meet the following criteria.
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
    mate status 0. Flags, tags and ranges may be specified in the
    \code{ScanBamParam} for fine tuning of results.
  } 

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

  \item{inter.feature}{
    A logical indicating if the counting \code{mode} should be aware of
    overlapping features. When TRUE (default), reads mapping to multiple
    features are dropped (i.e., not counted). When FALSE, these reads are
    retained and a count is assigned to each feature they map to.

    There are 6 possible combinations of the \code{mode} and
    \code{inter.feature} arguments. When \code{inter.feature=FALSE} the
    behavior of modes \sQuote{Union} and \sQuote{IntersectionStrict} are 
    essentially \sQuote{countOverlaps} with \sQuote{type=any} and
    \code{type=within}, respectively. \sQuote{IntersectionNotEmpty} does
    not reduce to a simple countOverlaps because common (shared) regions 
    of the annotation are removed before counting.}

  \item{singleEnd}{
    A logical value indicating if reads are single or paired-end.}

  \item{fragments}{
    A logical value indicating if singletons, reads with unmapped pairs
    and other fragments should be included in the counting. When 
    \code{fragments=FALSE} only reads paired with the algorithm described
    at ?findMateAlignment are counted. When \code{fragments=TRUE} (default)
    all singletons, reads with unmapped pairs and other fragments are
    counted in addition to the reads paired with the ?findMateAlignment 
    algorithm. This argument applies to paired-end reads only so
    \code{singleEnd} must be FALSE.}

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
    \item{asMates}{A logical(0) indicating if the records should be returned as
      mated pairs.}

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

    \item{isOpen}{Tests whether the \code{BamFile} \code{con} has been
      opened for reading.}

    \item{isIncomplete}{Tests whether the \code{BamFile} \code{con} is
      niether closed nor at the end of the file.}

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

    \item{asMates, asMates<-}{Return or set a logical(0)
      indicating if the records should be returned as mated pairs.}
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
  \code{\link[GenomicRanges:summarizeOverlaps]{summarizeOverlaps}} 
}

\examples{

##
## Iterating over Bam files.
##

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

##
## summarizeOverlaps().
##

## See ?'summarizeOverlaps' in the GenomicRanges package for examples.

##
## findSpliceOverlaps. 
##

## See ?'findSpliceOverlaps' for examples 
}

\keyword{classes}
