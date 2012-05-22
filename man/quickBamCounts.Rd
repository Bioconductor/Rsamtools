\name{quickBamCounts}

\alias{quickBamCounts}

\title{Group the records of a BAM file based on their flag bits and count
       the number of records in each group}

\description{
\code{quickBamCounts} groups the records of a BAM file based on their flag
bits and counts the number of records in each group.
}

\usage{
quickBamCounts(file, index=file, param=NULL, main.groups.only=FALSE)
}

\arguments{
  \item{file, index}{The path to the BAM file to read, and to the index
    file of the BAM file to read, respectively. The latter is given
    \emph{without} the '.bai' extension. See \code{\link{scanBam}} for more
    information.}

  \item{param}{\code{NULL} or an instance of \code{\linkS4class{ScanBamParam}}.
    This determines which records are considered in the counting.}

  \item{main.groups.only}{If \code{TRUE}, then the counting is performed
    for the main groups only.
  }
}

\value{
  Nothing is returned. A summary of the counts is printed to the console
  unless redirected by \code{\link{sink}}.
}

\author{H. Pages}

\references{
  \url{http://samtools.sourceforge.net/}
}

\seealso{
  \code{\link{scanBam}},
  \code{\link{ScanBamParam}}
}

\examples{
bamfile <- system.file("extdata", "ex1.bam", package="Rsamtools")
quickBamCounts(bamfile)
}

\keyword{manip}
