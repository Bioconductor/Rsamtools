\name{quickBamFlagSummary}

\alias{quickBamFlagSummary}
\alias{quickBamFlagSummary,character-method}
\alias{quickBamFlagSummary,list-method}

\title{Group the records of a BAM file based on their flag bits and count
       the number of records in each group}

\description{
\code{quickBamFlagSummary} groups the records of a BAM file based on their flag
bits and counts the number of records in each group.
}

\usage{
quickBamFlagSummary(file, ..., param=ScanBamParam(), main.groups.only=FALSE)

\S4method{quickBamFlagSummary}{character}(file, index=file, ..., param=ScanBamParam(),
    main.groups.only=FALSE)

\S4method{quickBamFlagSummary}{list}(file, ..., param=ScanBamParam(), main.groups.only=FALSE)
}

\arguments{

  \item{file, index}{For the character method, the path to the BAM file
    to read, and to the index file of the BAM file to read,
    respectively.

    For the list() method, \code{file} is a named list with elements
    \dQuote{qname} and \dQuote{flag} with content as from
    \code{\link{scanBam}}.}

  \item{...}{Additional arguments, perhaps used by methods.} 

  \item{param}{An instance of \code{\linkS4class{ScanBamParam}}.
    This determines which records are considered in the counting.}

  \item{main.groups.only}{If \code{TRUE}, then the counting is performed
    for the main groups only.
  }
}

\value{
  Nothing is returned. A summary of the counts is printed to the console
  unless redirected by \code{\link{sink}}.
}

\author{Hervé Pagès}

\references{
  \url{http://samtools.sourceforge.net/}
}

\seealso{
  \code{\link{scanBam}},
  \code{\link{ScanBamParam}}.

  \code{\link{BamFile}} for a method for that class.
 }

\examples{
bamfile <- system.file("extdata", "ex1.bam", package="Rsamtools",
                       mustWork=TRUE)
quickBamFlagSummary(bamfile)
}

\keyword{manip}
