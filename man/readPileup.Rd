\name{readPileup}
\Rdversion{1.1}
\alias{readPileup}
\alias{readPileup,connection-method}
\alias{readPileup,character-method}

\title{
  Import samtools 'pileup' files.
}
\description{
  Import files created by evaluation of samtools' \code{pileup -cv}
  command.
}
\usage{
readPileup(file, ...)
\S4method{readPileup}{connection}(file, ..., variant=c("SNP", "indel", "all"))
}

\arguments{

  \item{file}{The file name, or
	\code{\link[base:connections]{connection}}, of the pileup output file
	to be parsed.}

  \item{\dots}{Additional arguments, passed to methods. For instance,
	specify \code{variant} for the readPileup,character-method.}

  \item{variant}{Type of variant to parse; select one.}

}

\value{

  \code{readPileup} returns a \code{\link[GenomicRanges]{GRanges}} object.
  
  The value returned by \code{variant="SNP"} or \code{variant="all"}
  contains:
  \describe{

	\item{space:}{The chromosome names (fastq ids) of the reference
	  sequence}

	\item{position:}{The nucleotide position (base 1) of the variant.}

	\item{referenceBase:}{The nucleotide in the reference sequence.}

	\item{consensusBase;}{The consensus nucleotide, as determined by
	  samtools pileup.}

	\item{consensusQuality:}{The phred-scaled consensus quality.}

	\item{snpQuality:}{The phred-scaled SNP quality (probability of the
	  consensus being identical to the reference).}

	\item{maxMappingQuality:}{The root mean square mapping quality of reads
	  overlapping the site.}

	\item{coverage:}{The number of reads covering the site.}

  }

  The value returned by \code{variant="indel"} contains space, position,
  reference, consensus, consensusQuality, snpQuality, maxMappingQuality,
  and coverage fields, and:

  \describe{

	\item{alleleOne, alleleTwo}{The first (typically, in the reference
	  sequence) and second allelic variants.}

	\item{alleleOneSupport, alleleTwoSupport}{The number of reads
	  supporting each allele.}

	\item{additionalIndels}{The number of additional indels present.}

  }

}
\references{
  \url{http://samtools.sourceforge.net/}
}
\author{
  Sean Davis
}

\examples{

fl <- system.file("extdata", "pileup.txt", package="Rsamtools")
(res <- readPileup(fl))
xtabs(~referenceBase + consensusBase, mcols(res))[DNA_BASES,]

\dontrun{## uses a pipe, and arguments passed to read.table
## three successive piles of 100 records each
cmd <- "samtools pileup -cvf human_b36_female.fa.gz na19240_3M.bam"
p <- pipe(cmd, "r")
snp <- readPileup(p, nrow=100)  # variant="SNP"
indel <- readPileup(p, nrow=100, variant="indel")
all <- readPileup(p, nrow=100, variant="all")
}

}

\keyword{ manip }
