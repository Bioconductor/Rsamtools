\name{applyPileups}
\Rdversion{1.1}
\alias{applyPileups}
% \alias{yieldBam}

\title{

  Create summary pile-up statistics across multiple BAM files.

}
\description{

  \code{applyPileups} scans one or more BAM files, returning
  position-specific sequence and quality summaries.

}
\usage{
applyPileups(files, FUN, ..., param)
}

\arguments{

  \item{files}{A \code{\link{PileupFiles}} instances.}

  \item{FUN}{A function of 1 argument, \code{x}, to be evaluated for
    each yield (see \code{yieldSize}, \code{yieldBy},
    \code{yieldAll}). The argument \code{x} is a \code{list}, with
    elements describing the current pile-up. The elements of the list
    are determined by the argument \code{what}, and include:

    \describe{

      \item{seqnames:}{(Always returned) A named \code{integer()}
	representing the seqnames corresponding to each position
	reported in the pile-up.  This is a run-length encoding, where
	the names of the elements represent the seqnames, and the values
	the number of successive positions corresponding to that
	seqname.}

      \item{pos:}{Always returned) A \code{integer()} representing the
	genomic coordinate of each pile-up position.}

      \item{seq:}{An \code{array} of dimensions nucleotide x file x
	position.

	The \sQuote{nucleotide} dimension is length 5, corresponding to
	\sQuote{A}, \sQuote{C}, \sQuote{G}, \sQuote{T}, and \sQuote{N}
	respectively.

	Entries in the array represent the number of times the
	nucleotide occurred in reads in the file overlapping the
	position.

      }

      \item{qual:}{Like \code{seq}, but summarizing quality; the first
	dimension is the Phred-encoded quality score, ranging from
	\sQuote{!} (0) to \sQuote{~} (93).}

    }
  }

  \item{\dots}{Additional arguments, passed to methods.}

  \item{param}{An instance of the object returned by
    \code{PileupParam}.}

}

\value{

  \code{applyPileups} returns a \code{list} equal in length to the
  number of times \code{FUN} has been called, with each element
  containing the result of \code{FUN}.

  \code{PileupParam} returns an object describing the parameters.

}

\references{\url{http://samtools.sourceforge.net/}}

\author{Martin Morgan}

\seealso{\code{\link{PileupParam}}.}

\examples{

fl <- system.file("extdata", "ex1.bam", package="Rsamtools")

fls <- PileupFiles(c(fl, fl))

calcInfo <-
    function(x)
{
    ## information at each pile-up position
    info <- apply(x[["seq"]], 2, function(y) {
        y <- y[c("A", "C", "G", "T"),]
        y <- y + 1L                     # continuity
        cvg <- colSums(y)
        p <- y / cvg[col(y)]
        h <- -colSums(p * log(p))
        ifelse(cvg == 4L, NA, h)
    })
    list(seqnames=x[["seqnames"]], pos=x[["pos"]], info=info)
}
which <- GRanges(c("seq1", "seq2"), IRanges(c(1000, 1000), 2000))
param <- PileupParam(which=which, what="seq")
res <- applyPileups(fls, calcInfo, param=param)
str(res)
head(res[[1]][["pos"]])		# positions matching param
head(res[[1]][["info"]])		# inforamtion in each file

## 'param' as part of 'files'
fls1 <- PileupFiles(c(fl, fl), param=param)
res1 <- applyPileups(fls1, calcInfo)
identical(res, res1)

## yield by position, across ranges
param <- PileupParam(which=which, yieldSize=500L, yieldBy="position",
                     what="seq")
res <- applyPileups(fls, calcInfo, param=param)
sapply(res, "[[", "seqnames")

}

\keyword{ manip }
