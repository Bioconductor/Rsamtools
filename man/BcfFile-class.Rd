\name{BcfFile}
\Rdversion{1.1}
\docType{class}
\alias{BcfFile-class}
\alias{BcfFileList-class}
% con/destructors
\alias{BcfFile}
\alias{BcfFileList}
\alias{open.BcfFile}
\alias{close.BcfFile }
% accessors
\alias{isOpen,BcfFile-method}
\alias{bcfMode}
% methods
\alias{scanBcfHeader,BcfFile-method}
\alias{scanBcf,BcfFile-method}
\alias{indexBcf,BcfFile-method}

\title{Manipulate BCF files.}

\description{

  Use \code{BcfFile()} to create a reference to a BCF (and optionally
  its index). The reference remains open across calls to methods, 
  avoiding costly index re-loading.

  \code{BcfFileList()} provides a convenient way of managing a list of
  \code{BcfFile} instances.

}

\usage{

## Constructors

BcfFile(file, index = file,
        mode=ifelse(grepl("\\\\.bcf$", file), "rb", "r"))
BcfFileList(...)

## Opening / closing

\S3method{open}{BcfFile}(con, ...)
\S3method{close}{BcfFile}(con, ...)

## accessors; also path(), index()

\S4method{isOpen}{BcfFile}(con, rw="")
bcfMode(object)

## actions

\S4method{scanBcfHeader}{BcfFile}(file, ...)
\S4method{scanBcf}{BcfFile}(file, ..., param=ScanBcfParam())
\S4method{indexBcf}{BcfFile}(file, ...)

}

\arguments{

  \item{con, object}{An instance of \code{BcfFile}.}

  \item{file}{A character(1) vector of the BCF file path or, (for
    indexBcf) an instance of \code{BcfFile} point to a BCF file.}

  \item{index}{A character(1) vector of the BCF index.}

  \item{mode}{A character(1) vector; \code{mode="rb"} indicates a binary
    (BCF) file, \code{mode="r"} a text (VCF) file.}

  \item{param}{An optional \code{\linkS4class{ScanBcfParam}} instance to
     further influence scanning.}

  \item{...}{Additional arguments. For \code{BcfFileList}, this can
    either be a single character vector of paths to BCF files, or
    several instances of \code{BcfFile} objects.}

  \item{rw}{Mode of file; ignored.}

}

\section{Objects from the Class}{

  Objects are created by calls of the form \code{BcfFile()}.

}

\section{Fields}{

  The \code{BcfFile} class inherits fields from the
  \code{\linkS4class{RsamtoolsFile}} class.

}

\section{Functions and methods}{

  \code{BcfFileList} inherits methods from
  \code{\link{RsamtoolsFileList}} and \code{\link{SimpleList}}.

  Opening / closing:
  \describe{

    \item{open.BcfFile}{Opens the (local or remote) \code{path} and
      \code{index} (if \code{bamIndex} is not \code{character(0)}),
      files.  Returns a \code{BcfFile} instance.}

    \item{close.BcfFile}{Closes the \code{BcfFile} \code{con}; returning
	  (invisibly) the updated \code{BcfFile}. The instance may be
	  re-opened with \code{open.BcfFile}.}

  }

  Accessors: 
  \describe{

    \item{path}{Returns a character(1) vector of the BCF path name.}

    \item{index}{Returns a character(1) vector of BCF index name.}

    \item{bcfMode}{Returns a character(1) vector BCF mode.}

  }

  Methods:
  \describe{

    \item{scanBcf}{Visit the path in \code{path(file)}, returning the
      result of \code{\link{scanBcf}} applied to the specified path.}

    \item{show}{Compactly display the object.}

  }

}
\author{Martin Morgan}

\examples{

fl <- system.file("extdata", "ex1.bcf", package="Rsamtools")
bf <- BcfFile(fl)        # implicit index
bf
identical(scanBcf(bf), scanBcf(fl))

rng <- GRanges(c("seq1", "seq2"), IRanges(1, c(1575, 1584)))
param <- ScanBcfParam(which=rng)
bcf <- scanBcf(bf, param=param)  ## all ranges

## ranges one at a time 'bf'
open(bf)
sapply(seq_len(length(rng)), function(i, bcfFile, rng) {
    param <- ScanBcfParam(which=rng)
    bcf <- scanBcf(bcfFile, param=param)[[1]]
    ## do extensive work with bcf
    isOpen(bf)  ## file remains open
}, bf, rng)

}

\keyword{classes}
