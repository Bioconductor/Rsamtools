\name{BamViews}
\Rdversion{1.1}
\docType{class}
\alias{BamViews-class}
\alias{BamViews}
% accessors
\alias{bamPaths}
\alias{bamSamples}
\alias{bamRanges}
\alias{bamExperiment}
% methods
\alias{dim,BamViews-method}
\alias{[,BamViews,ANY,missing-method}
\alias{[,BamViews,missing,ANY-method}
\alias{[,BamViews,ANY,ANY-method}
\alias{[,BamViews,ANY,missing-method}
\alias{[,BamViews,RangesList,ANY-method}
\alias{[,BamViews,RangesList,missing-method}
\alias{show,BamViews-method}

\title{Views into a set of BAM files}
\description{

  Use \code{BamViews()} to reference a set of disk-based BAM files to be
  processed (e.g., queried using \code{\link{scanBam}}) as a single
  \sQuote{experiment}.

}

\usage{

## Constructor
BamViews(bamPaths=character(0),
         bamSamples=new("DataFrame", nrows=length(bamPaths)),
         bamRanges = RangedData(), bamExperiment = list(),
         ...)
## Accessors
bamPaths(object)
bamSamples(object)
bamRanges(object)
bamExperiment(object)

## Subset
\S4method{[}{BamViews,ANY,ANY}(x, i, j, \dots, drop=TRUE)
\S4method{[}{BamViews,ANY,missing}(x, i, j, \dots, drop=TRUE)
\S4method{[}{BamViews,missing,ANY}(x, i, j, \dots, drop=TRUE)

## Show
\S4method{show}{BamViews}(object)

}

\arguments{

  \item{bamPaths}{A character() vector of path names.}

  \item{bamSamples}{A \code{\linkS4class{DataFrame}} instance with as
	many rows as \code{length(bamPaths)}, containing sample information
	associated with each path.}

  \item{bamRanges}{A \code{\linkS4class{RangedData}} instance with
	ranges defined on the spaces of the BAM files. Ranges are \emph{not}
	validated against the BAM files.}

  \item{bamExperiment}{A list() containing additional information about
	the experiment.}

  \item{...}{Additional arguments, currently unused.}

  \item{x}{An instance of \code{BamViews}.}

  \item{object}{An instance of \code{BamViews}.}

  \item{i}{During subsetting, a logical or numeric index into
	\code{bamRanges}.}

  \item{j}{During subsetting, a logical or numeric index into
	\code{bamSamples} and \code{bamPaths}.}

  \item{drop}{A logical(1), \emph{ignored} by all \code{BamViews}
	subsetting methods.}

}

\section{Objects from the Class}{

  Objects are created by calls of the form \code{BamViews()}.

}
\section{Slots}{
  \describe{

  \item{bamPaths}{A character() vector of path names.}

  \item{bamSamples}{A \code{\linkS4class{DataFrame}} instance with as
	many rows as \code{length(bamPaths)}, containing sample information
	associated with each path.}

  \item{bamRanges}{A \code{\linkS4class{RangedData}} instance with
	ranges defined on the spaces of the BAM files. Ranges are \emph{not}
	validated against the BAM files.}

  \item{bamExperiment}{A list() containing additional information about
	the experiment.}

  }
}

\section{Functions and methods}{

  See 'Usage' for details on invocation.

  Constructor:
  \describe{

	\item{BamViews:}{Returns a \code{BamViews} object.}

  }

  Accessors:
  \describe{

	\item{bamPaths}{Returns a character() vector of path names.}

	\item{bamSamples}{Returns a \code{\linkS4class{DataFrame}} instance
	  with as many rows as \code{length(bamPaths)}, containing sample
	  information associated with each path.}

	\item{bamRanges}{Returns a \code{\linkS4class{RangedData}} instance
	  with ranges defined on the spaces of the BAM files. Ranges are
	  \emph{not} validated against the BAM files.}

	\item{bamExperiment}{Returns a list() containing additional
	  information about the experiment.}

  }

  Methods:
  \describe{

	\item{"["}{Subset the object by \code{bamRanges} or
	  \code{bamSamples}.}

	\item{show}{Compactly display the object.}

  }

}
\author{Martin Morgan}

\examples{

fls <- list.files(system.file("extdata", package="Rsamtools"),
                  "\\\\.bam$", full=TRUE)
rl <- RangesList(chr1=IRanges(seq(10000, 90000, 10000), width=500),
                 chr2=IRanges(seq(100000, 900000, 100000), width=5000))
v <- BamViews(fls, bamRanges=RangedData(rl, Count=seq_len(18L)))
v
v[1:5,]
bamRanges(v[c(1:5, 11:15),])
}
\keyword{classes}
