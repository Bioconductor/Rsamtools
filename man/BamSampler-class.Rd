\name{BamSampler}
\Rdversion{1.1}
\docType{class}
\alias{BamSampler-class}
% con/destructors
\alias{BamSampler}
% methods
\alias{scanBam,BamSampler-method}
\alias{show,BamSampler-method}

\title{Sample from a BAM files}

\description{

  Use \code{BamSampler()} to create a reference to a BAM file (and
  optionally its index). Calls to \code{scanBam} (and many functions
  that use \code{scanBam}) draw a random sample from the BAM file.

}

\usage{

## Constructors

BamSampler(file, index = file, ..., yieldSize, obeyQname = FALSE, asMates = FALSE) 

\S4method{scanBam}{BamSampler}(file, index=file, ..., param=ScanBamParam(what=scanBamWhat()))
}

\arguments{

  \item{file}{character(1); BAM file path for \code{BamSampler}, or
  \code{BamSampler} index for \code{scanBam} and other functions.}

  \item{index}{character(1); the BAM index file path (for
    \code{BamFile}); ignored for other methods.}

  \item{...}{Additional arguments; see \code{\linkS4class{BamFile}-class}.}

  \item{yieldSize}{integer(1); number of records to yield each time the
    file is read from using \code{scanBam}.}

  \item{obeyQname}{logical(1); indicating whether the file is sorted by
    \code{qname} and if so, that qnames are not split between yields.}

  \item{asMates}{logical(1); indicating whether the records should be
    returned as mated pairs.}

  \item{param}{An optional \code{\linkS4class{ScanBamParam}} instance to
     further influence scanning, counting, or filtering.}

}

\section{Objects from the Class}{

  Objects are created by calls of the form \code{BamSampler()}.

}

\section{Fields}{

  The \code{BamSampler} class inherits fields from the
  \code{\linkS4class{BamFile}} class.

}

\section{Functions and methods}{

  \code{BamSampler} inherits methods from \code{\link{BamFile}} and can
  be used in place of \code{BamFile} in many functions.

}

\author{Martin Morgan}

\examples{
library(GenomicAlignments)  # for readGAlignmentsFromBam()
fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
samp <- BamSampler(fl, yieldSize=1000)
## two independent samples
head(readGAlignmentsFromBam(samp))
head(readGAlignmentsFromBam(samp))

}

\keyword{classes}
