\name{pileupBam}
\Rdversion{1.1}
\alias{.PileupParam}
\alias{.pileupBam}

\title{

  Create summary pile-up statistics across multiple BAM files.

}
\description{

  \code{pileupBam} scans one or more BAM files, returning
  position-specific sequence and quality summaries.

  \code{PileupParam} creates parameters to influence pile-up
  attributes.

}
\usage{
.pileupBam(files, FUN=identity, ..., param)
.PileupParam(flag=scanBamFlag(), minBaseQuality=13L,
    minMapQuality=0L, minDepth=0L, maxDepth=250L, yieldSize=1L,
    yieldBy=c("range", "position"), yieldAll=FALSE,
    which=GRanges(), what=c("seq", "qual"))
}

\arguments{

  \item{files}{A \code{list()} of open \code{\link{BamFile}} instances.}

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

	The \sQuote{nucleotide} dimension is length 16,
	corresponding to the samtools nucleotide encoding, in particular
	elements 2, 3, 5, and 9 correspond to nucleotides \sQuote{A},
	\sQuote{C}, \sQuote{G}, and \sQuote{T}, respectively.

	Entries in the array represent the number of times the
	nucleotide occurred in reads in the file overlapping the
	position.

      }

      \item{qual:}{Like \code{seq}, but summarizing quality; the first
	dimension is the Phred-encoded quality score.}

    }
  }

  \item{\dots}{Additional arguments, passed to methods.}

  \item{param}{An instance of the object returned by
    \code{PileupParam}.}

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

}

\value{

  \code{pileupBam} returns a \code{list} equal in length to the number
  of times \code{FUN} has been called, with each element containing the
  result of \code{FUN}.

  \code{PileupParam} returns an object describing the parameters.

}
\references{
  \url{http://samtools.sourceforge.net/}
}
\author{
  Martin Morgan
}

\examples{

fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
PileupParam <- Rsamtools:::.PileupParam
pileupBam <- Rsamtools:::.pileupBam

fls <- list(open(BamFile(fl)), open(BamFile(fl)))

calcInfo <-
    function(x)
{
    ## information at each pile-up position
    x[[3]] <- apply(x[[3]], 2, function(y) {
        y <- y[c(1, 2, 4, 8) + 1,]      # A, C, G, T
        y <- y + 1L                     # continuity
        cvg <- colSums(y)
        p <- y / cvg[col(y)]
        h <- -colSums(p * log(p))
        ifelse(cvg == 4L, NA, h)
    })
    x
}
which <- GRanges(c("seq1", "seq2"), IRanges(c(1000, 1000), 2000))
param <- PileupParam(which=which)
res <- pileupBam(fls, calcInfo, param=param)
str(res)
head(res[[1]][["pos"]])		# positions matching param
head(res[[1]][["seq"]])		# inforamtion in each file

## yield by position, across ranges
param <- PileupParam(which=which, yieldSize=2000L, yieldBy="position")
res <- pileupBam(fls, calcInfo, param=param)
str(res)
res[[1]][["seqnames"]]

}

\keyword{ manip }
