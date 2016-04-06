\name{BamViews}
\Rdversion{1.1}
\docType{class}
\alias{BamViews-class}
\alias{BamViews}
\alias{BamViews,missing-method}
\alias{BamViews,GRanges-method}
% accessors
\alias{bamPaths}
\alias{bamIndicies}
\alias{bamSamples}
\alias{bamRanges}
\alias{bamExperiment}
\alias{bamDirname<-}
\alias{bamSamples<-}
\alias{bamRanges<-}
% methods
\alias{names,BamViews-method}
\alias{names<-,BamViews-method}
\alias{dim,BamViews-method}
\alias{dimnames,BamViews-method}
\alias{dimnames<-,BamViews-method}
\alias{dimnames<-,BamViews,ANY-method}
\alias{[,BamViews,ANY,missing-method}
\alias{[,BamViews,missing,ANY-method}
\alias{[,BamViews,ANY,ANY-method}
\alias{[,BamViews,ANY,missing-method}
\alias{scanBam,BamViews-method}
\alias{countBam,BamViews-method}
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
     bamIndicies=bamPaths,
     bamSamples=DataFrame(row.names=make.unique(basename(bamPaths))),
     bamRanges, bamExperiment = list(), ...)
\S4method{BamViews}{missing}(bamPaths=character(0),
     bamIndicies=bamPaths,
     bamSamples=DataFrame(row.names=make.unique(basename(bamPaths))),
     bamRanges, bamExperiment = list(), ..., auto.range=FALSE)
## Accessors
bamPaths(x)
bamSamples(x)
bamSamples(x) <- value
bamRanges(x)
bamRanges(x) <- value
bamExperiment(x)

\S4method{names}{BamViews}(x)
\S4method{names}{BamViews}(x) <- value
\S4method{dimnames}{BamViews}(x)
\S4method{dimnames}{BamViews,ANY}(x) <- value

bamDirname(x, ...) <- value

## Subset
\S4method{[}{BamViews,ANY,ANY}(x, i, j, \dots, drop=TRUE)
\S4method{[}{BamViews,ANY,missing}(x, i, j, \dots, drop=TRUE)
\S4method{[}{BamViews,missing,ANY}(x, i, j, \dots, drop=TRUE)

## Input
\S4method{scanBam}{BamViews}(file, index = file, ..., param = ScanBamParam(what=scanBamWhat()))
\S4method{countBam}{BamViews}(file, index = file, ..., param = ScanBamParam())

## Show
\S4method{show}{BamViews}(object)
}

\arguments{

  \item{bamPaths}{A character() vector of BAM path names.}

  \item{bamIndicies}{A character() vector of BAM index file path names,
    \emph{without} the \sQuote{.bai} extension.}

  \item{bamSamples}{A \code{\linkS4class{DataFrame}} instance with as
	many rows as \code{length(bamPaths)}, containing sample information
	associated with each path.}

  \item{bamRanges}{Missing or a \code{\linkS4class{GRanges}} instance with
        ranges defined on the reference chromosomes of the BAM files. Ranges
        are \emph{not} validated against the BAM files.}

  \item{bamExperiment}{A list() containing additional information about
	the experiment.}

  \item{auto.range}{If \code{TRUE} and all \code{bamPaths} exist,
    populate the ranges with the union of ranges returned in the
    \code{target} element of \code{scanBamHeader}.}

  \item{...}{Additional arguments.}

  \item{x}{An instance of \code{BamViews}.}

  \item{object}{An instance of \code{BamViews}.}

  \item{value}{An object of appropriate type to replace content.}

  \item{i}{During subsetting, a logical or numeric index into
	\code{bamRanges}.}

  \item{j}{During subsetting, a logical or numeric index into
	\code{bamSamples} and \code{bamPaths}.}

  \item{drop}{A logical(1), \emph{ignored} by all \code{BamViews}
	subsetting methods.}

  %% input

  \item{file}{An instance of \code{BamViews}.}

  \item{index}{A character vector of indices, corresponding to the
     \code{bamPaths(file)}. }

  \item{param}{An optional \code{\linkS4class{ScanBamParam}} instance to
     further influence scanning or counting.}
}

\section{Objects from the Class}{

  Objects are created by calls of the form \code{BamViews()}.

}
\section{Slots}{
  \describe{

  \item{bamPaths}{A character() vector of BAM path names.}

  \item{bamIndicies}{A character() vector of BAM index path names.}

  \item{bamSamples}{A \code{\linkS4class{DataFrame}} instance with as
	many rows as \code{length(bamPaths)}, containing sample information
	associated with each path.}

  \item{bamRanges}{A \code{\linkS4class{GRanges}} instance with
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

	\item{bamPaths}{Returns a character() vector of BAM path names.}

	\item{bamIndicies}{Returns a character() vector of BAM index path
      names.}

	\item{bamSamples}{Returns a \code{\linkS4class{DataFrame}} instance
	  with as many rows as \code{length(bamPaths)}, containing sample
	  information associated with each path.}

	\item{bamSamples<-}{Assign a \code{\linkS4class{DataFrame}} instance
	  with as many rows as \code{length(bamPaths)}, containing sample
	  information associated with each path.}

	\item{bamRanges}{Returns a \code{\linkS4class{GRanges}} instance
	  with ranges defined on the spaces of the BAM files. Ranges are
	  \emph{not} validated against the BAM files.}

	\item{bamRanges<-}{Assign a \code{\linkS4class{GRanges}} instance
	  with ranges defined on the spaces of the BAM files. Ranges are
	  \emph{not} validated against the BAM files.}

	\item{bamExperiment}{Returns a list() containing additional
	  information about the experiment.}

    \item{names}{Return the column names of the \code{BamViews}
      instance; same as \code{names(bamSamples(x))}.}

    \item{names<-}{Assign the column names of the \code{BamViews}
      instance.}

    \item{dimnames}{Return the row and column names of the
      \code{BamViews} instance.}

    \item{dimnames<-}{Assign the row and column names of the
      \code{BamViews} instance.}

  }

  Methods:
  \describe{

    \item{"["}{Subset the object by \code{bamRanges} or \code{bamSamples}.}

    \item{scanBam}{Visit each path in \code{bamPaths(file)}, returning
      the result of \code{scanBam} applied to the specified
      path. \code{bamRanges(file)} takes precedence over
      \code{bamWhich(param)}.}

    \item{countBam}{Visit each path in \code{bamPaths(file)}, returning
      the result of \code{countBam} applied to the specified
      path. \code{bamRanges(file)} takes precedence over
      \code{bamWhich(param)}.}

    \item{show}{Compactly display the object.}

  }

}
\author{Martin Morgan}

\examples{
fls <- system.file("extdata", "ex1.bam", package="Rsamtools",
                   mustWork=TRUE)
rngs <- GRanges(seqnames = Rle(c("chr1", "chr2"), c(9, 9)),
                ranges = c(IRanges(seq(10000, 90000, 10000), width=500),
                           IRanges(seq(100000, 900000, 100000), width=5000)),
                Count = seq_len(18L))
v <- BamViews(fls, bamRanges=rngs)
v
v[1:5,]
bamRanges(v[c(1:5, 11:15),])
bamDirname(v) <- getwd()
v
}

\keyword{classes}
