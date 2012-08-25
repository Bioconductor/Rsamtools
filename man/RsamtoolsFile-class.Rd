\name{RsamtoolsFile}
\Rdversion{1.1}
\docType{class}
\alias{RsamtoolsFile-class}
\alias{index}
\alias{path}
\alias{yieldSize}
\alias{yieldSize<-}
\alias{yieldSize<-,RsamtoolsFile-method}
\alias{path,RsamtoolsFile-method}
\alias{isOpen,RsamtoolsFile-method}
\alias{show,RsamtoolsFile-method}
\alias{yieldSize,RsamtoolsFile-method}

\title{A base class for managing file references in Rsamtools}

\description{

  \code{RsamtoolsFile} is a base class for managing file references in
  \pkg{Rsamtools}; it is not intended for direct use by users -- see, e.g.,
  \code{\link{BamFile}}.

}

\usage{

## accessors
index(object)
\S4method{path}{RsamtoolsFile}(object, ...)
\S4method{isOpen}{RsamtoolsFile}(con, rw="")
\S4method{yieldSize}{RsamtoolsFile}(object, ...)
yieldSize(object, ...) <- value
\S4method{show}{RsamtoolsFile}(object)

}

\arguments{

  \item{con, object}{An instance of a class derived from \code{RsamtoolsFile}.}

  \item{rw}{Mode of file; ignored.}

  \item{\dots}{Additional arguments, unused.}

  \item{value}{Replacement value.}
}

\section{Objects from the Class}{

  Users do not directly create instances of this class; see, e.g.,
  \code{\link{BamFile}-class}.

}

\section{Fields}{

  The \code{RsamtoolsFile} class is implemented as an S4 reference
  class. It has the following fields:

  \describe{

    \item{.extptr}{An \code{externalptr} initialized to an internal
      structure with opened bam file and bam index pointers.}

    \item{path}{A character(1) vector of the file name.}

    \item{index}{A character(1) vector of the index file name.}

    \item{yieldSize}{An integer(1) vector of the number of records to
      yield.} 

  }
}

\section{Functions and methods}{

  Accessors:
  \describe{

    \item{path}{Returns a character(1) vector of path names.}

    \item{index}{Returns a character(1) vector of index path names.}

    \item{yieldSize, yieldSize<-}{Return or set an integer(1) vector
      indicating yield size.}

  }

  Methods:
  \describe{

    \item{isOpen}{Report whether the file is currently open.}

    \item{show}{Compactly display the object.}

  }

}
\author{Martin Morgan}

\keyword{classes}
