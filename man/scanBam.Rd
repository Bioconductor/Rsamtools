\name{BamInput}
\Rdversion{1.1}

\alias{scanBam}
\alias{scanBam,character-method}
\alias{countBam}
\alias{countBam,character-method}
\alias{scanBamHeader}
\alias{scanBamHeader,character-method}
\alias{indexBam}
\alias{indexBam,character-method}
\alias{filterBam}
\alias{filterBam,character-method}

\title{

  Import, count, index, and other operations on `BAM' (binary alignment)
  files.

}
\description{
  Import binary `BAM' files into a list structure, with facilities for
  selecting what fields and which records are imported.
}
\usage{

scanBam(file, index=file, ..., param=ScanBamParam())

countBam(file, index=file, ..., param=ScanBamParam())

scanBamHeader(files, ...)
\S4method{scanBamHeader}{character}(files, ...)

indexBam(files, ...)
\S4method{indexBam}{character}(files, ...)

filterBam(file, destination, index=file, ...)
\S4method{filterBam}{character}(file, destination, index=file, ...,
    indexDestination=TRUE, param=ScanBamParam())

}

\arguments{

  \item{file}{The character(1) file name of the `BAM' file to be processed.}

  \item{files}{The character() file names of the `BAM' file to be indexed.}

  \item{index}{The character(1) name of the index file of the 'BAM' file
	being processed; this is given \emph{without} the '.bai' extension.}

  \item{destination}{The character(1) file name of the location where
    the filtered output file will be created.}

  \item{\dots}{Additional arguments, passed to methods.}

  \item{indexDestination}{A logical(1) indicating whether the created
    destination file should also be indexed.}

  \item{param}{An instance of \code{\linkS4class{ScanBamParam}}. This
	influences what fields and which records are imported.}

}
\details{

  The \code{scanBam} function parses binary BAM files; text SAM files
  can be parsed using R's \code{\link{scan}} function, especially with
  arguments \code{what} to control the fields that are parsed.

  The \code{countBam} function returns a count of records consistent
  with \code{param}.

  The \code{scanBamHeader} function visits the header information in a
  BAM file, returning for each file a list containing elements
  \code{targets} and \code{text}, as described below. The SAM / BAM
  specification does not require that the content of the header be
  consistent with the content of the file, e.g., more targets may be
  present that are represented by reads in the file.

  The \code{indexBam} function creates an index for each BAM file
  specified, analogous to the \sQuote{samtools index} function.

  The \code{filterBam} function parses records in \code{file} satisfying
  the \code{bamWhich} of \code{param}, writing each record satisfying
  the \code{bamFlag} and \code{bamSimpleCigar} criteria of \code{param}
  to file \code{destination}. An index file is created on the
  destination when \code{indexDestination=TRUE}.
  
  Details of the \code{ScanBamParam} class are provide on its help page;
  several salient points are reiterated here. \code{ScanBamParam} can
  contain a field \code{what}, specifying the components of the BAM
  records to be returned. Valid values of \code{what} are available with
  \code{\link{scanBamWhat}}. \code{ScanBamParam} can contain an argument
  \code{which} that specifies a subset of reads to return. This requires
  that the BAM file be indexed, and that the file be named following
  samtools convention as \code{<bam_filename>.bai}. \code{ScanBamParam}
  can contain an argument \code{tag} to specify which tags will be
  extracted. 

}

\value{

  The \code{scanBam,character-method} returns a list of lists. The outer
  list groups results from each \code{Ranges} list of
  \code{bamWhich(param)}; the outer list is of length one when
  \code{bamWhich(param)} has length 0. Each inner list contains elements
  named after \code{scanBamWhat()}; elements omitted from
  \code{bamWhat(param)} are removed. The content of non-null elements
  are as follows, taken from the description in the samtools API
  documentation:

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

  \item{qwidth}{The width of the query, as calculated from the
	\code{cigar} encoding; normally equal to the width of the query
	returned in \code{seq}.}

  \item{mrnm}{The reference to which the mate (of a paired end or mate
	pair read) aligns.}

  \item{mpos}{The position to which the mate aligns.}

  \item{isize}{Inferred insert size for paired end alignments.}

  \item{seq}{The query sequence, in the 5' to 3' orientation. If aligned
	to the minus strand, it is the reverse complement of the original
	sequence.}

  \item{qual}{Phred-encoded, phred-scaled base quality score, oriented
	as \code{seq}.}

  \code{scanBamHeader} returns a list, with one element for each file
  named in \code{files}. The list contains two element. The
  \code{targets} element contains target (reference) sequence
  lengths. The \code{text} element is itself a list with each element a
  list corresponding to tags (e.g., \sQuote{@SQ}) found in the header,
  and the associated tag values.

  \code{indexBam} returns the file name of the index file created.

  \code{filterBam} returns the file name of the destination file
  created.

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
table(width(res0[["seq"]])) # query widths
table(res0[["qwidth"]], useNA="always") # query widths derived from cigar
table(res0[["cigar"]], useNA="always")
table(res0[["strand"]], useNA="always")
table(res0[["flag"]], useNA="always")

which <- RangesList(seq1=IRanges(1000, 2000),
                    seq2=IRanges(c(100, 1000), c(1000, 2000)))
p1 <- ScanBamParam(which=which)
res1 <- scanBam(fl, param=p1)
names(res1)
names(res1[[2]])

p2 <- ScanBamParam(what=c("rname", "strand", "pos", "qwidth"))
res2 <- scanBam(fl, param=p2)
                
p3 <- ScanBamParam(flag=scanBamFlag(isMinusStrand=FALSE))
length(scanBam(fl, param=p3)[[1]])

}

\keyword{ manip }
