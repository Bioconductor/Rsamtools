\name{scanBam}
\Rdversion{1.1}
\alias{scanBam}
\alias{scanBam,character-method}

\title{
  Import samtools `BAM' (binary) files.
}
\description{
  Import binary `BAM' files into a list structure, with facilities for
  selecting which fields and what records are to be imported.
}
\usage{
scanBam(file, index=file, ...)
\S4method{scanBam}{character}(file, index=file, ..., param=ScanBamParam())
}

\arguments{

  \item{file}{The file name of the `BAM' file to be parsed.}

  \item{index}{The name of the index file of the 'BAM' file being
	parsed; this is given \emph{without} the '.bai' extension.}

  \item{\dots}{Additional arguments, passed to methods.}

  \item{param}{For the \code{scanBam,character-method}, an instance of
	\code{\linkS4class{ScanBamParam}}. This influences which fields and
	what records are input.}

}
\details{

  This method parses binary BAM files; text SAM files can be parsed
  using R's \code{\link{scan}} function, especially with arguments
  \code{what} to control the fields are parsed.
  
  Details of the \code{ScanBamParam} class are provide on its help page;
  several salient points are reiterated here. \code{ScanBamParam} can
  contain a field \code{what}, specifying the components of the BAM
  records to be returned. Valid values of \code{what} are available with
  \code{\link{scanBamWhat}}. \code{ScanBamParam} can contain an
  argument \code{which} that specifies a subset of reads to return. This
  requires that the BAM file be indexed, and that the file be named
  following samtools convention as \code{<bam_filename>.bai}.

}

\value{

  The \code{scanBam,character-method} returns a list of lists. The outer
  list groups results from each \code{Ranges} list of
  \code{bamWhich(param)}; the outer list is of length one when
  \code{bamWhich(param)} has length 0. Each inner list contains elements
  named after \code{scanBamWhat()}; elements omitted from
  \code{bamWhat(param)} have value \code{NULL}. The content of non-null
  elements are as follows, taken from the description in the samtools
  API documentation:

  \item{qname}{The query name, i.e., identifier, associated with the
	read.}

  \item{flag}{A numeric value summarizing details of the read. See
	\code{\linkS4class{ScanBamParam}} and the \code{flag} argument, and
	\code{scanBamFlag()}.}

  \item{rname}{The name of the reference to which the read is aligned.}

  \item{strand}{The strand to which the read is aligned.}

  \item{pos}{The genomic coordinate at the start of the
	alignment. Coordinates are `left-most', i.e., at the 3' end of a
	read on the '-' strand, and 1-based. The position \emph{excludes}
	clipped nucleotides, even though soft-clipped nucleotides are
	included in \code{seq}.}

  \item{width}{The width of the read, as calculated from the
	\code{cigar} encoding; normally equal to the width of the aligned
	read returned in \code{seq}.}

  \item{mrnm}{The reference to which the mate (of a paired end
	read) aligns.}

  \item{mpos}{The position to which the mate aligns.}

  \item{isize}{Inferred insert size for paired end alignments.}

  \item{seq}{The read sequence, in its original (unaligned, 5' to 3')
	orientation.}

  \item{qual}{Phred-encoded, phred-scaled base quality score, oriented
	as \code{seq}.}

}
\references{
  \url{http://samtools.sourceforge.net/}
}
\author{
  Martin Morgan
}

\seealso{

  \code{\linkS4class{ScanBamParam}}, \code{\link{scanBamWhat}},
  \code{\link{scanBamFlag}}

}
\examples{

fl <- system.file("extdata", "ex1.bam", package="Rsamtools")

res0 <- scanBam(fl)[[1]] # always list-of-lists
names(res0)
length(res0[["qname"]])
lapply(res0, head, 3)
table(width(res0[["seq"]])) # read widths
table(res0[["width"]], useNA="always") # aligned width
table(cigars(res0[["cigar"]]), useNA="always")
table(res0[["strand"]], useNA="always")
table(res0[["flag"]], useNA="always")

which <- RangesList(seq1=IRanges(1000, 2000),
                    seq2=IRanges(c(100, 1000), c(1000, 2000)))
p1 <- ScanBamParam(which=which)
res1 <- scanBam(fl, param=p1)
names(res1)
names(res1[[2]])

p2 <- ScanBamParam(what=c("rname", "strand", "pos", "width"))
res2 <- scanBam(fl, param=p2)
                
p3 <- ScanBamParam(flag=scanBamFlag(isMinusStrand=FALSE))
length(scanBam(fl, param=p3)[[1]])

}

\keyword{ manip }
