\name{BcfFile}
\Rdversion{1.1}
\docType{class}
\alias{BcfFile-class}
% con/destructors
\alias{BcfFile}
\alias{open.BcfFile}
\alias{close.BcfFile }
% accessors
\alias{isOpen,BcfFile-method}
\alias{bcfMode}
% methods
\alias{scanBcfHeader,BcfFile-method}
\alias{scanBcf,BcfFile-method}

\title{Manipulate BCF or VCF files.}

\description{

  Use \code{BcfFile()} to create a reference to a BCF (and optionally
  its index) or VCF file. The reference remains open across calls to
  methods, avoiding costly index re-loading.

}

\usage{

## Opening / closing

BcfFile(file, index = file,
        mode=ifelse(grepl("\\\\.bcf$", file), "rb", "r"))
\S3method{open}{BcfFile}(con, ...)
\S3method{close}{BcfFile}(con, ...)

## accessors; also path(), index()

\S4method{isOpen}{BcfFile}(con, rw="")
bcfMode(object)

## actions

\S4method{scanBcfHeader}{BcfFile}(file, ...)
\S4method{scanBcf}{BcfFile}(file, ..., param=ScanBcfParam())

}

\arguments{

  \item{con, object}{An instance of \code{BcfFile}.}

  \item{file}{A character(1) vector of the VCF or BCF file path.}

  \item{index}{A character(1) vector of the BCF index.}

  \item{mode}{A character(1) vector; \code{mode="rw"} indicates a binary
    (BCF) file, \code{mode="r"} a text (VCF) file.}

  \item{param}{An optional \code{\linkS4class{ScanBcfParam}} instance to
     further influence scanning.}

  \item{...}{Additional arguments, currently unused.}

  \item{rw}{Mode of file; ignored.}

}

\section{Objects from the Class}{

  Objects are created by calls of the form \code{BcfFile()}.

}

\section{Fields}{

  The \code{BamFile} class inherits fields from the
  \code{\linkS4class{RsamtoolsFile}} class.

}

\section{Functions and methods}{

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
