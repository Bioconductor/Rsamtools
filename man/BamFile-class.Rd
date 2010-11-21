\name{BamFile}
\Rdversion{1.1}
\docType{class}
\alias{BamFile-class}
% con/destructors
\alias{openBam}
\alias{openBam,character-method}
\alias{openBam,BamFile-method}
\alias{close.BamFile }
\alias{show,BamFile-method}
% accessors
\alias{bamPath}
\alias{bamIndex}
\alias{isOpen,BamFile-method}
% methods
\alias{scanBamHeader,BamFile-method}
\alias{scanBam,BamFile-method}
\alias{countBam,BamFile-method}
\alias{filterBam,BamFile-method}
\alias{indexBam,BamFile-method}
\alias{sortBam,BamFile-method}
\alias{readBamGappedAlignments,BamFile-method}

\title{Maintain open BAM files across calls to methods}

\description{

  Use \code{openBam()} to create a reference to a BAM file (and
  optionally its index). The reference can avoid costly index loading
  across calls.

}

\usage{

## Opening / closing / showing
openBam(bamPath, bamIndex=bamPath, ...)
\S4method{openBam}{character}(bamPath, bamIndex=bamPath, ...)
\S4method{openBam}{BamFile}(bamPath, bamIndex=bamPath, ...)
\S3method{close}{BamFile}(con, ...)
\S4method{show}{BamFile}(object)

## accessors
bamPath(object)
bamIndex(object)
\S4method{isOpen}{BamFile}(con, rw="")

## actions
\S4method{scanBamHeader}{BamFile}(files, ...)
\S4method{scanBam}{BamFile}(file, index=file, ..., param=ScanBamParam())
\S4method{countBam}{BamFile}(file, index=file, ..., param=ScanBamParam())
\S4method{filterBam}{BamFile}(file, destination, index=file, ...,
    indexDestination=TRUE, param=ScanBamParam())
\S4method{indexBam}{BamFile}(files, ...)
\S4method{sortBam}{BamFile}(file, destination, ..., byQname=FALSE, maxMemory=512)
\S4method{readBamGappedAlignments}{BamFile}(file, index, ..., which)

}

\arguments{

  \item{...}{Additional arguments.}

  \item{con, object, file, files}{An instance of \code{BamFile}.}

  \item{index}{A character vector of indices; ignored for all methods on
    this page. }

  \item{bamPath}{The file name (local or remote) of the BAM file.}

  \item{bamIndex}{The file name (local or remote) of the BAM index
    file.}

  \item{destination}{character(1) file path to write filtered reads to.}

  \item{indexDestination}{logical(1) indicating whether the destination
	file should also be indexed.}

  \item{byQname, maxMemory}{See \code{\link{sortBam}}.}

  \item{param}{An optional \code{\linkS4class{ScanBamParam}} instance to
     further influence scanning, counting, or filtering.}

  \item{which}{An optional \code{\linkS4class{RangesList}} instance to
     further subset \code{file}.}

  \item{rw}{Mode of file; ignored.}

}

\section{Objects from the Class}{

  Objects are created by calls of the form \code{openBam()}.

}

\section{Fields}{

  The \code{BamFile} class is implemented as a S4-reference class. It
  has the following fields:

  \describe{

    \item{.extptr}{An \code{externalptr} initialized to an internal
      structure with opened bam file and bam index pointers.}

    \item{bamPath}{A character(1) vector of the BAM file name.}

    \item{bamIndex}{A character(1) vector of the BAM index file name.}

  }
}

\section{Functions and methods}{

  See 'Usage' for details on invocation.

  Opening / closing:
  \describe{

    \item{openBam}{Opens the (local or remote) \code{bamPath} and
      \code{bamIndex} (if \code{bamIndex} is not \code{character(0)}),
      files.  Returns a \code{BamFile} object.}

	\item{closeBam}{Close the \code{BamFile} \code{con}; returning
	  (invisibly) the updated \code{BamFile}. The object may be
	  re-opened with \code{openBam}.}

  }

  Accessors:
  \describe{

    \item{bamPath}{Returns a character(1) vector of BAM path names.}

    \item{bamIndex}{Returns a character(1) vector of BAM index path
      names.}

  }

  Methods:
  \describe{

    \item{scanBam}{Visit the path in \code{bamPath(file)}, returning the
      result of \code{\link{scanBam}} applied to the specified path.}

    \item{countBam}{Visit the path in \code{bamPath(file)}, returning
      the result of \code{\link{countBam}} applied to the specified
      path.}

    \item{filterBam}{Visit the path in \code{bamPath(file)}, returning
      the result of \code{\link{filterBam}} applied to the specified
      path.}

    \item{indexBam}{Visit the path in \code{bamPath(file)}, returning
      the result of \code{\link{indexBam}} applied to the specified
      path.}

    \item{sortBam}{Visit the path in \code{bamPath(file)}, returning the
      result of \code{\link{sortBam}} applied to the specified path.}

    \item{readBamGappedAlignments}{Visit the path in
      \code{bamPath(file)}, returning the result of
      \code{readBamGappedAlignments} applied to the specified path. See
      \code{\link{readBamGappedAlignments}}.}

    \item{show}{Compactly display the object.}

  }

}
\author{Martin Morgan}

\examples{

fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
bf <- openBam(fl)        # implicit bamIndex
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
