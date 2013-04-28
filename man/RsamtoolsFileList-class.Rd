\name{RsamtoolsFileList}
\Rdversion{1.1}
\docType{class}
\alias{RsamtoolsFileList-class}
\alias{RsamtoolsFileList}
\alias{path,RsamtoolsFileList-method}
\alias{isOpen,RsamtoolsFileList-method}
\alias{open.RsamtoolsFileList}
\alias{close.RsamtoolsFileList}
\alias{names,RsamtoolsFileList-method}
\alias{yieldSize,RsamtoolsFileList-method}
\alias{yieldSize<-,RsamtoolsFileList-method}

\title{A base class for managing lists of Rsamtools file references}

\description{

  \code{RsamtoolsFileList} is a base class for managing lists of file
  references in \pkg{Rsamtools}; it is not intended for direct use --
  see, e.g., \code{\link{BamFileList}}.

}

\usage{

\S4method{path}{RsamtoolsFileList}(object, ...)
\S4method{isOpen}{RsamtoolsFileList}(con, rw="")
\S3method{open}{RsamtoolsFileList}(con, ...)
\S3method{close}{RsamtoolsFileList}(con, ...)
\S4method{names}{RsamtoolsFileList}(x)
\S4method{yieldSize}{RsamtoolsFileList}(object, ...)

}

\arguments{

  \item{con, object, x}{An instance of a class derived from
    \code{RsamtoolsFileList}.}

  \item{rw}{Mode of file; ignored.}

  \item{\dots}{Additional arguments.}

}

\section{Objects from the Class}{

  Users do not directly create instances of this class; see, e.g.,
  \code{\link{BamFileList}-class}.

}

\section{Functions and methods}{

  This class inherits functions and methods for subseting, updating, and
  display from the \code{\link{SimpleList}} class.

  Methods:
  \describe{

    \item{isOpen:}{Report whether each file in the list is currently
      open.}

    \item{open:}{Attempt to open each file in the list.}

    \item{close:}{Attempt to close each file in the list.}

    \item{names:}{Names of each element of the list or, if names are
      \code{NULL}, the basename of the path of each element.}

  }

}
\author{Martin Morgan}

\keyword{classes}
