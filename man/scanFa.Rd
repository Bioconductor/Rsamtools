\name{FaInput}
\Rdversion{1.1}

\alias{indexFa}
\alias{indexFa,character-method}
\alias{scanFaIndex}
\alias{scanFaIndex,character-method}
\alias{countFa}
\alias{countFa,character-method}
\alias{scanFa}
\alias{scanFa,character,GRanges-method}
\alias{scanFa,character,RangesList-method}
\alias{scanFa,character,RangedData-method}
\alias{scanFa,character,missing-method}

\title{

  Operations on indexed 'fasta' files.

}
\description{

  Scan indexed fasta (or compressed fasta) files and their indicies.

}
\usage{

indexFa(file, ...)
\S4method{indexFa}{character}(file, ...)

scanFaIndex(file, ...)
\S4method{scanFaIndex}{character}(file, ...)

countFa(file, ...)
\S4method{countFa}{character}(file, ...)

scanFa(file, param, ...)
\S4method{scanFa}{character,GRanges}(file, param, ...)
\S4method{scanFa}{character,RangesList}(file, param, ...)
\S4method{scanFa}{character,RangedData}(file, param, ...)
\S4method{scanFa}{character,missing}(file, param, ...)

}

\arguments{

  \item{file}{A character(1) vector containing the fasta file path.}

  \item{param}{An optional \code{\linkS4class{GRanges}},
    \code{\linkS4class{RangesList}}, or \code{\linkS4class{RangedData}}
    instance to select reads (and sub-sequences) for input.}

  \item{...}{Additional arguments, currently unused.}

}

\value{

  \code{indexFa} visits the path in \code{file} and create an index file
  at the same location but with extension \sQuote{.fai}).

  \code{scanFaIndex} reads the sequence names and and widths of recorded
  in an indexed fasta file, returning the information as a
  \code{\linkS4class{GRanges}} object.

  \code{countFa} returns the number of records in the fasta file.

  \code{scanFa} return the sequences indicated by \code{param} as a
  \code{\linkS4class{DNAStringSet}} instance. \code{seqnames(param)}
  selects the sequences to return; \code{start(param)} and
  \code{end{param}} define the (1-based) region of the sequence to
  return. Values of \code{end(param)} greater than the width of the
  sequence are set to the width of the sequence. When \code{param} is
  missing, all records are selected. When \code{param} is
  \code{GRanges()}, no records are selected.

}

\references{

  \url{http://samtools.sourceforge.net/} provides information on
  \code{samtools}.

}

\author{

  Martin Morgan <mtmorgan@fhcrc.org>.

}

\examples{
fa <- system.file("extdata", "ce2dict1.fa", package="Rsamtools")
countFa(fa)
(idx <- scanFaIndex(fa))
(dna <- scanFa(fa, idx[1:2]))
ranges(idx) <- narrow(ranges(idx), -10)  # last 10 nucleotides
(dna <- scanFa(fa, idx[1:2]))
}

\keyword{ manip }
