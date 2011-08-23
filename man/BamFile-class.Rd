\name{BamFile}
\Rdversion{1.1}
\docType{class}
\alias{BamFile-class}
\alias{BamFileList-class}
% con/destructors
\alias{BamFile}
\alias{BamFileList}
\alias{open.BamFile}
\alias{close.BamFile }
% accessors
\alias{isOpen,BamFile-method}
% methods
\alias{scanBamHeader,BamFile-method}
\alias{seqinfo,BamFile-method}
\alias{scanBam,BamFile-method}
\alias{countBam,BamFile-method}
\alias{countBam,BamFileList-method}
\alias{filterBam,BamFile-method}
\alias{indexBam,BamFile-method}
\alias{sortBam,BamFile-method}
\alias{readBamGappedAlignments,BamFile-method}
\alias{readBamGappedReads,BamFile-method}

\title{Maintain SAM and BAM files}

\description{

  Use \code{BamFile()} to create a reference to a BAM file (and
  optionally its index). The reference remains open across calls to
  methods, avoiding costly index re-loading.

  \code{BamFileList()} provides a convenient way of managing a list of
  \code{BamFile} instances.

}

\usage{

## Constructors

BamFile(file, index=file)
BamFileList(...)

## Opening / closing

\S3method{open}{BamFile}(con, ...)
\S3method{close}{BamFile}(con, ...)

## accessors; also path(), index()

\S4method{isOpen}{BamFile}(con, rw="")

## actions

\S4method{scanBamHeader}{BamFile}(files, ...)
\S4method{seqinfo}{BamFile}(x)
\S4method{scanBam}{BamFile}(file, index=file, ..., param=ScanBamParam())
\S4method{countBam}{BamFile}(file, index=file, ..., param=ScanBamParam())
\S4method{countBam}{BamFileList}(file, index=file, ..., param=ScanBamParam())
\S4method{filterBam}{BamFile}(file, destination, index=file, ...,
    indexDestination=TRUE, param=ScanBamParam())
\S4method{indexBam}{BamFile}(files, ...)
\S4method{sortBam}{BamFile}(file, destination, ..., byQname=FALSE, maxMemory=512)
\S4method{readBamGappedAlignments}{BamFile}(file, index=file, use.names=FALSE, param=NULL)
\S4method{readBamGappedReads}{BamFile}(file, index=file, use.names=FALSE, param=NULL)

}

\arguments{

  \item{...}{Additional arguments. For \code{BamFileList}, this can
    either be a single character vector of paths to BAM files, or
    several instances of \code{BamFile} objects.}

  \item{con}{An instance of \code{BamFile}.}

  \item{x, file, files}{A character vector of BAM file paths (for
    \code{BamFile}) or a \code{BamFile} instance (for other methods).}

  \item{index}{A character vector of indices (for \code{BamFile});
    ignored for all other methods on this page.}

  \item{destination}{character(1) file path to write filtered reads to.}

  \item{indexDestination}{logical(1) indicating whether the destination
	file should also be indexed.}

  \item{byQname, maxMemory}{See \code{\link{sortBam}}.}

  \item{param}{An optional \code{\linkS4class{ScanBamParam}} instance to
     further influence scanning, counting, or filtering.}

  \item{use.names}{Construct the names of the returned object
     from the query template names (QNAME field)? If not (the default),
     then the returned object has no names.}

  \item{which}{An optional \code{\linkS4class{RangesList}} instance to
     further subset \code{file}.}

  \item{rw}{Mode of file; ignored.}

}

\section{Objects from the Class}{

  Objects are created by calls of the form \code{BamFile()}.

}

\section{Fields}{

  The \code{BamFile} class inherits fields from the
  \code{\linkS4class{RsamtoolsFile}} class.

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

    \item{readBamGappedAlignments, readBamGappedReads}{Visit
      the path in \code{path(file)}, returning the result of
      \code{readBamGappedAlignments} or \code{readBamGappedReads}
      applied to the specified path.
      See \code{\link{readBamGappedAlignments}}.}

    \item{show}{Compactly display the object.}

  }

}
\author{Martin Morgan}

\examples{

fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
bf <- open(BamFile(fl))        # implicit index
bf
identical(scanBam(bf), scanBam(fl))

rng <- GRanges(c("seq1", "seq2"), IRanges(1, c(1575, 1584)))

## repeatedly visit 'bf'
sapply(seq_len(length(rng)), function(i, bamFile, rng) {
    param <- ScanBamParam(which=rng[i], what="seq")
    bam <- scanBam(bamFile, param=param)[[1]]
    alphabetFrequency(bam[["seq"]], baseOnly=TRUE, collapse=TRUE)
}, bf, rng)

}

\keyword{classes}
