\name{yieldReduce}
\alias{yieldReduce}

\title{
  Iterate through a BAM (or other) file, reducing output to a single result.
}
\description{

  Rsamtools files can be created with a \sQuote{yieldSize} argument that
  influences the number of records (chunk size) input at one time (see,
  e.g,. \code{\link{BamFile}}). \code{yieldReduce} iterates through the
  file, processing each chunk and reducing it with previously input
  chunks. This is a memory efficient way to process large data files,
  especially when the final result fits in memory.

}
\usage{
yieldReduce(X, MAP, REDUCE, DONE, ..., init, ITERATE = TRUE)
}

\arguments{

  \item{X}{A \code{\link{BamFile}} instance (or other class for which
    \code{isOpen}, \code{open}, \code{close} methods are defined, and
    which support input of sequential chunks).}

  \item{MAP}{A function of one or more arguments, \code{X, ...},
    returning a \code{VALUE} passed to \code{DONE} and \code{REDUCE}.}

  \item{REDUCE}{A function of one (\code{ITERATE=FALSE} or two
    (\code{ITERATE=TRUE}) arguments, returning the reduction (e.g.,
    addition) of the argument(s). If missing, \code{REDUCE} is \code{c}
    (when \code{ITERATE=TRUE}) or \code{identity} when
    (when \code{ITERATE=FALSE}).}

  \item{DONE}{A function of one argument, the \code{VALUE} of the most
    recent call to \code{MAP(X, ...)}. If missing, \code{DONE} is
    \code{function(VALUE) length(VALUE) == 0}.
  }

  \item{\dots}{Additional arguments, passed to \code{MAP}.}

  \item{init}{(Optional) Initial value used for \code{REDUCE} when
    \code{ITERATE=TRUE}.}

  \item{ITERATE}{logical(1) determining whether the call to
    \code{REDUCE} is iterative (\code{ITERATE=TRUE}) or cumulative
    (\code{ITERATE=FALSE}).

  }
}
\details{
  
  When \code{ITERATE=TRUE}, \code{REDUCE} is initially invoked with
  either the \code{init} value and the value of the first call to
  \code{MAP} or, if \code{init} is missing, the values of the first two
  calls to \code{MAP}.

  When \code{ITERATE=FALSE}, \code{REDUCE} is invoked with a list
  containing a list with as many elements as there were calls to
  \code{MAP}.  Each element the result of an invocation of \code{MAP}.

}
\value{

  The return value is the value returned by the final invocation of
  \code{REDUCE}, or \code{init} if provided and no data were yield'ed,
  or \code{list()} if \code{init} is missing and no data were yield'ed.

}

\author{Martin Morgan \url{mtmorgan@fhcrc.org}}

\seealso{
  \code{\link{BamFile}}, \code{\link{TabixFile}}, \code{\linkS4class{RsamtoolsFile}}.
}
\examples{
fl <- system.file(package="Rsamtools", "extdata", "ex1.bam")

## nucleotide frequency of mapped reads
bf <- BamFile(fl, yieldSize=500) ## typically, yieldSize=1e6
param <- ScanBamParam(
    flag=scanBamFlag(isUnmappedQuery=FALSE),
    what="seq")
MAP <- function(X, param) {
    value <- scanBam(X, param=param)[[1]][["seq"]]
    if (length(value))
        alphabetFrequency(value, collapse=TRUE)
    else value       # will be integer(0)
}
REDUCE <- `+`        # add successive alphabetFrequency matrices 
yieldReduce(bf, MAP, REDUCE, param=param)

## coverage
if (require(GenomicAlignments)) {
    MAP <- function(X)
        coverage(readGAlignments(X))
    REDUCE <- `+`
    DONE <- function(VALUE)
        ## coverage() on zero GAlignments returns an RleList,
        ## each element of  which has 0 coverage
        sum(sum(VALUE)) == 0L
    yieldReduce(bf, MAP, REDUCE, DONE)
}

}
\keyword{manip}
