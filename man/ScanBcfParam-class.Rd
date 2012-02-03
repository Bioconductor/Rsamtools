\name{ScanBcfParam-class}
\Rdversion{1.1}
\docType{class}

\alias{ScanBVcfParam-class}

\alias{ScanBcfParam}
\alias{ScanBcfParam-class}

\alias{ScanBcfParam,GRanges-method}
\alias{ScanBcfParam,RangedData-method}
\alias{ScanBcfParam,RangesList-method}
\alias{ScanBcfParam,missing-method}
\alias{bcfInfo}
\alias{bcfGeno}
\alias{bcfTrimEmpty}
\alias{bcfWhich}

\alias{show,ScanBVcfParam-method}

\title{Parameters for scanning BCF files}
\description{

  Use \code{ScanBcfParam()} to create a parameter object influencing the
  \sQuote{INFO} and \sQuote{GENO} fields parsed, and which recrods are
  imported from a BCF file. Use of \code{which} requires that a BCF
  index file (\code{<filename>.bci}) exists.

}

\usage{

ScanBcfParam(info=character(), geno=character(), trimEmpty=TRUE,
             which, asGRanges=FALSE, ...)
\S4method{ScanBcfParam}{missing}(info=character(), geno=character(), trimEmpty=TRUE,
             which, asGRanges=FALSE, ...)
\S4method{ScanBcfParam}{RangesList}(info=character(), geno=character(), trimEmpty=TRUE,
             which, asGRanges=FALSE, ...)
\S4method{ScanBcfParam}{RangedData}(info=character(), geno=character(), trimEmpty=TRUE,
             which, asGRanges=FALSE, ...)
\S4method{ScanBcfParam}{GRanges}(info=character(), geno=character(), trimEmpty=TRUE,
             which, asGRanges=FALSE, ...)

## Accessors
bcfInfo(object)
bcfGeno(object)
bcfTrimEmpty(object)
bcfWhich(object)

}

\arguments{

  \item{info}{A character() vector of \sQuote{INFO} fields (see
    \link[VariantAnnotation]{scanVcfHeader}) to be returned.}

  \item{geno}{A character() vector of \sQuote{GENO} fields (see
    \link[VariantAnnotation]{scanVcfHeader}) to be returned. \code{character(0)}
    returns all fields, \code{NA_character_} returns none.}

  \item{trimEmpty}{A logical(1) indicating whether \sQuote{GENO} fields
    with no values should be returned.}

  \item{which}{An object, for which a method is defined (see usage,
    above), describing the sequences and ranges to be queried. Variants
    whose \code{POS} lies in the interval(s) \code{[start, end)} are
    returned.} 

  \item{asGRanges}{A logical(1) for use with \code{ScanVcfParam} only.}

  \item{object}{An instance of class \code{ScanBcfParam}.}

  \item{\dots}{Arguments used internally.}

}

\section{Objects from the Class}{

  Objects can be created by calls of the form \code{ScanBcfParam()}.

}
\section{Slots}{
  \describe{
    \item{\code{which}:}{Object of class \code{"RangesList"} indicating
      which reference sequence and coordinate variants must overlap.}
    \item{\code{info}:}{Object of class \code{"character"} indicating
      portions of \sQuote{INFO} to be returned.}
    \item{\code{geno}:}{Object of class \code{"character"} indicating
      portions of \sQuote{GENO} to be returned. }
    \item{\code{trimEmpty}:}{Object of class \code{"logical"} indicating
      whether empty \sQuote{GENO} fields are to be returned.}
  }
}
\section{Functions and methods}{

  See 'Usage' for details on invocation.

  Constructor:
  \describe{

    \item{ScanBcfParam:}{Returns a \code{ScanBcfParam} object. The 
      \code{which} argument to the constructor can be one of several types, 
      as documented above.}  }

  Accessors:
  \describe{

    \item{bcfInfo, bcfGeno, bcfTrimEmpty, bcfWhich:}{Return the
      corresponding field from \code{object}.}

  }

  Methods:
  \describe{

    \item{show}{Compactly display the object.}

  }
}
\author{
  Martin Morgan \url{mtmorgan@fhcrc.org}
}
\seealso{

  \link[VariantAnnotation]{scanVcf}

}

\examples{
## see ?ScanVcfParam examples
}

\keyword{classes}
