\name{PileupFiles}
\Rdversion{1.1}
\docType{class}
\alias{PileupFiles-class}
% con/destructors
\alias{PileupFiles}
\alias{PileupFiles,character-method}
\alias{PileupFiles,list-method}
\alias{open.PileupFiles}
\alias{close.PileupFiles}
% accessors
\alias{isOpen,PileupFiles-method}
\alias{plpFiles}
\alias{plpParam}
% methods
\alias{applyPileups,PileupFiles,PileupParam-method}
\alias{applyPileups,PileupFiles,missing-method}
\alias{show,PileupFiles-method}

\title{Represent BAM files for pileup summaries.}

\description{

  Use \code{PileupFiles()} to create a reference to a BAM files (and
  their indicies), to be used for calculating pile-up summaries.

}

\usage{

## Constructors
PileupFiles(files, ..., param=PileupParam())
\S4method{PileupFiles}{character}(files, ..., param=PileupParam())
\S4method{PileupFiles}{list}(files, ..., param=PileupParam())

## opening / closing
\S3method{open}{PileupFiles}(con, ...)
\S3method{close}{PileupFiles}(con, ...)

## accessors; also path()
\S4method{isOpen}{PileupFiles}(con, rw="")
plpFiles(object)
plpParam(object)

## actions
\S4method{applyPileups}{PileupFiles,missing}(files, FUN, ..., param)
\S4method{applyPileups}{PileupFiles,PileupParam}(files, FUN, ..., param)

## display
\S4method{show}{PileupFiles}(object)

}

\arguments{

  \item{files}{For \code{PileupFiles}, a \code{character()} or
    \code{list} of \code{BamFile} instances representing files to be
    included in the pileup. Using a \code{list} of \code{BamFile} allows
    indicies to be specified when these are in non-standard format. All
    elements of \code{\ldots} must be the same type.

    For \code{applyPileups,PileupFiles-method}, a \code{PileupFiles}
    instance.}

  \item{\ldots}{Additional arguments, currently ignored.}

  \item{con, object}{An instance of \code{PileupFiles}.}

  \item{FUN}{A function of one argument; see \code{\link{applyPileups}}.}

  \item{param}{An instance of \code{\link{PileupParam}}, 
    to select which records to include in the pileup, and which summary
    information to return.}

  \item{rw}{character() indicating mode of file; not used for
    \code{TabixFile}.}

}

\section{Objects from the Class}{

  Objects are created by calls of the form \code{PileupFiles()}.

}

\section{Fields}{

  The \code{PileupFiles} class is implemented as an S4 reference
  class. It has the following fields:

  \describe{

    \item{files}{A list of \code{\link{BamFile}} instances.}

    \item{param}{An instance of \code{\link{PileupParam}}.}

  }
}

\section{Functions and methods}{

  Opening / closing:
  \describe{

    \item{open.PileupFiles}{Opens the (local or remote) \code{path} and
      \code{index} of each file in the \code{PileupFiles}
      instance. Returns a \code{PileupFiles} instance.}

    \item{close.PileupFiles}{Closes each file in the \code{PileupFiles}
      instance; returning (invisibly) the updated
      \code{PileupFiles}. The instance may be re-opened with
      \code{open.PileupFiles}.}

  }

  Accessors: 
  \describe{

    \item{plpFiles}{Returns the \code{list} of the files in the
      \code{PileupFiles} instance.}

    \item{plpParam}{Returns the \code{\link{PileupParam}} content of the
      \code{PileupFiles} instance.}

  }

  Methods:
  \describe{

    \item{applyPileups}{Calculate the pileup across all files in
      \code{files} according to criteria in \code{param} (or
      \code{plpParam(files)} if \code{param} is missing), invoking
      \code{FUN} on each range or collection of positions. See
      \code{\link{applyPileups}}.}

    \item{show}{Compactly display the object.}

  }

}

\author{Martin Morgan}

\examples{
example(applyPileups)
}

\keyword{classes}
