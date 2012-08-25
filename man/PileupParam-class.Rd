\name{PileupParam}
\Rdversion{1.1}
\docType{class}
\alias{PileupParam-class}
\alias{PileupParam}
% accessors
\alias{plpFlag<-}
\alias{plpFlag}
\alias{plpMaxDepth<-}
\alias{plpMaxDepth}
\alias{plpMinBaseQuality<-}
\alias{plpMinBaseQuality}
\alias{plpMinDepth<-}
\alias{plpMinDepth}
\alias{plpMinMapQuality<-}
\alias{plpMinMapQuality}
\alias{plpWhat<-}
\alias{plpWhat}
\alias{plpWhich<-}
\alias{plpWhich}
\alias{plpYieldAll<-}
\alias{plpYieldAll}
\alias{plpYieldBy<-}
\alias{plpYieldBy}
\alias{plpYieldSize<-}
\alias{plpYieldSize}
% methods
\alias{show,PileupParam-method}

\title{Parameters for creating pileups from BAM files}
\description{

  Use \code{PileupParam()} to create a parameter object influencing
  what fields and which records  are used to calculate pile-ups, and to
  influence the values returned.

}

\usage{

# Constructor
PileupParam(flag = scanBamFlag(),
    minBaseQuality = 13L, minMapQuality = 0L,
    minDepth = 0L, maxDepth = 250L,
    yieldSize = 1L, yieldBy = c("range", "position"), yieldAll = FALSE,
    which = GRanges(), what = c("seq", "qual"))

# Accessors
plpFlag(object)
plpFlag(object) <- value
plpMaxDepth(object)
plpMaxDepth(object) <- value
plpMinBaseQuality(object)
plpMinBaseQuality(object) <- value
plpMinDepth(object)
plpMinDepth(object) <- value
plpMinMapQuality(object)
plpMinMapQuality(object) <- value
plpWhat(object)
plpWhat(object) <- value
plpWhich(object)
plpWhich(object) <- value
plpYieldAll(object)
plpYieldAll(object) <- value
plpYieldBy(object)
plpYieldBy(object) <- value
plpYieldSize(object)
plpYieldSize(object) <- value

\S4method{show}{PileupParam}(object)

}

\arguments{

  \item{flag}{An instance of the object returned by
    \code{\link{scanBamFlag}}, restricting various aspects of reads to
    be included or excluded.}

  \item{minBaseQuality}{The minimum read base quality below which the
    base is ignored when summarizing pileup information.}

  \item{minMapQuality}{The minimum mapping quality below which the
    entire read is ignored.}

  \item{minDepth}{The minimum depth of the pile-up below which the
    position is ignored.}

  \item{maxDepth}{The maximum depth of reads considered at any position;
    this can be used to limit memory consumption.}

  \item{yieldSize}{The number of records to include in each call to
    \code{FUN}.}

  \item{yieldBy}{How records are to be counted. By range (in which case
    \code{yieldSize} must equal 1) means that \code{FUN} is invoked once
    for each range in \code{which}. By position means that \code{FUN} is
    invoked whenever pile-ups have been accumulated for \code{yieldSize}
    positions, regardless of ranges in \code{which}.}

  \item{yieldAll}{Whether to report all positions
    (\code{yieldAll=TRUE}), or just those passing the filtering criteria
    of \code{flag}, \code{minBaseQuality}, etc. When
    \code{yieldAll=TRUE}, positions not passing filter criteria have
    \sQuote{0} entries in \code{seq} or \code{qual}.}

  \item{which}{A \code{GRanges} or \code{RangesList} instance
    restricting pileup calculations to the corresponding genomic
    locations.}

  \item{what}{A \code{character()} instance indicating what values are
    to be returned. One or more of \code{c("seq", "qual")}.}

  \item{object}{An instace of class \code{PileupParam}.}

  \item{value}{An instance to be assigned to the corresponding slot of
    the \code{PileupParam} instance.}

}

\section{Objects from the Class}{

  Objects are created by calls of the form \code{PileupParam()}.

}
\section{Slots}{

  Slot interpretation is as described in the \sQuote{Arguments} section.

  \describe{

    \item{\code{flag}}{Object of class \code{integer} encoding flags
      to be kept when they have their '0' (\code{keep0}) or '1'
      (\code{keep1}) bit set.}

    \item{\code{minBaseQuality}}{An \code{integer(1)}.}

    \item{\code{minMapQuality}}{An \code{integer(1)}.}

    \item{\code{minDepth}}{An \code{integer(1)}.}

    \item{\code{maxDepth}}{An \code{integer(1)}.}

    \item{\code{yieldSize}}{An \code{integer(1)}.}

    \item{\code{yieldBy}}{An \code{character(1)}.}

    \item{\code{yieldAll}}{A \code{logical(1)}.}

    \item{\code{which}}{A \code{GRanges} or \code{RangesList} instance.}

    \item{\code{what}}{A \code{character()}.}

  }
}

\section{Functions and methods}{

  See 'Usage' for details on invocation.

  Constructor:
  \describe{

    \item{PileupParam:}{Returns a \code{PileupParam} object.}

  }

  Accessors: get or set corresponding slot values; for setters,
  \code{value} is coerced to the type of the corresponding slot.

  \describe{

    \item{plpFlag, plpFlag<-}{Returns or sets the named \code{integer}
      vector of flags; see \code{\link{scanBamFlag}}.}

    \item{plpMinBaseQuality, plpMinBaseQuality<-}{Returns or sets an
      \code{integer(1)} vector of miminum base qualities.}

    \item{plpMinMapQuality, plpMinMapQuality<-}{Returns or sets an
      \code{integer(1)} vector of miminum map qualities.}

    \item{plpMinDepth, plpMinDepth<-}{Returns or sets an
      \code{integer(1)} vector of miminum pileup depth.}

    \item{plpMaxDepth, plpMaxDepth<-}{Returns or sets an
      \code{integer(1)} vector of the maximum depth to which pileups are
      calculated.}

    \item{plpYieldSize, plpYieldSize<-}{Returns or sets an
      \code{integer(1)} vector of yield size.}

    \item{plpYieldBy, plpYieldBy<-}{Returns or sets an
      \code{character(1)} vector determining how pileups will be
      returned.}

    \item{plpYieldAll, plpYieldAll<-}{Returns or sets an
      \code{logical(1)} vector indicating whether all positions, or just
      those satisfying pileup positions, are to be returned.}

    \item{plpWhich, plpWhich<-}{Returns or sets the object influencing
      which locations pileups are calculated over.}

    \item{plpWhat, plpWhat<-}{Returns or sets the \code{character}
      vector describing what summaries are returned by pileup.}

  }

  Methods:
  \describe{

	\item{show}{Compactly display the object.}

  }

}
\author{Martin Morgan}

\seealso{\link{applyPileups}.}

\examples{
example(applyPileups)
}
\keyword{classes}
