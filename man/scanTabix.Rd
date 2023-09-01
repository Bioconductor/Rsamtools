\name{TabixInput}
\Rdversion{1.1}

\alias{scanTabix}
\alias{scanTabix,character,GRanges-method}
\alias{scanTabix,character,IntegerRangesList-method}

\title{

  Operations on `tabix' (indexed, tab-delimited) files.

}
\description{
  Scan compressed, sorted, tabix-indexed, tab-delimited files.
}
\usage{

scanTabix(file, ..., param)
\S4method{scanTabix}{character,IntegerRangesList}(file, ..., param)
\S4method{scanTabix}{character,GRanges}(file, ..., param)

}

\arguments{

  \item{file}{The character() file name(s) of the tabix file be
    processed, or more flexibly an instance of class
    \code{\link{TabixFile}}.}

  \item{param}{A instance of \code{GRanges} or \code{IntegerRangesList}
    providing the sequence names and regions to be parsed.}


  \item{...}{Additional arguments, currently ignored.}

}

\value{

  \code{scanTabix} returns a list, with one element per region. Each element
  of the list is a character vector representing records in the region. If
  \code{param} is empty then all records will be returned.

}


\section{Error}{

  \code{scanTabix} signals errors using \code{signalCondition}. The
  following errors are signaled:

  \describe{

    \item{\code{scanTabix_param}}{\code{yieldSize(file)} must be NA when
      more than one range is specified.}

    \item{\code{scanTabix_io}}{A read error occured while inputing the
      tabix file. This might be because the file is corrupt, or of
      incorrect format (e.g., when \code{path} points to a plain text
      file but \code{index} is present, implying that \code{path} should
      be a \code{bgzip}ed file. The error message may include an error
      code representing the logical OR of these cryptic signals: 1,
      BGZF_ERR_ZLIB; 2, BGZF_ERR_HEADER; 4, BGZF_ERR_IO; 8,
      BGZF_ERR_MISUSE.}

  }
}

\references{

  \url{http://samtools.sourceforge.net/tabix.shtml}

}

\author{

  Martin Morgan <mtmorgan@fhcrc.org>.

}

\examples{
example(TabixFile)
}

\keyword{ manip }
