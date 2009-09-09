\name{ScanBamParam}
\Rdversion{1.1}
\docType{class}
\alias{ScanBamParam-class}
\alias{ScanBamParam}
% helpers
\alias{scanBamWhat}
\alias{scanBamFlag}
% accessors
\alias{bamWhat}
\alias{bamWhich}
\alias{bamFlag}
\alias{bamSimpleCigar}
% methods
\alias{show,ScanBamParam-method}

\title{Parameters for scanning BAM files}
\description{
  
  Use \code{ScanBamParam()} to create a parameter object influencing
  which reads and what fields are imported from a (binary) BAM file. Use
  of \code{which} requires that an BAM index file
  (\code{<filename>.bai}) exists.

}

\usage{

# Constructor
ScanBamParam(flag = scanBamFlag(), simpleCigar = FALSE, what =
    scanBamWhat(), which = RangesList())

# Constructor helpers
scanBamFlag (isPaired = NA, isProperPair = NA, isUnmappedQuery = NA, 
    hasUnmappedMate = NA, isMinusStrand = NA, isMateMinusStrand = NA, 
    isFirstMateRead = NA, isSecondMateRead = NA, isPrimaryRead = NA, 
    isValidVendorRead = NA, isDuplicate = NA) 
scanBamWhat()

# Accessors
bamWhat(object)
bamWhich(object)
bamFlag(object)
bamSimpleCigar(object)

\S4method{show}{ScanBamParam}(object)

}

\arguments{

  \item{flag}{An integer(2) vector used to filter reads based on their
	'flag' entry. This is most easily created with the
	\code{scamBamFlag()} helper function.}

  \item{simpleCigar}{A logical(1) vector which, when TRUE, returns only
	those reads for which the cigar (run-length encoded representation
	of the alignment) is missing or contains only matches / mismatches
	(\code{'M'}).}

  \item{what}{A character vector naming the fields to
	return. \code{scanBamWhat()} returns a vector of available
	fields. Fields are described on the \code{\link{scanBam}} help
	page.}

  \item{which}{A \code{\linkS4class{RangesList}} object, consisting of
	named \code{Ranges} elements. Names correspond to reference
	sequences, and ranges to the regions on that reference sequence for
	which matches are desired. Only reads overlapping the specified
	ranges are returned.}

  \item{isPaired}{A logical(1) indicating whether unpaired (FALSE),
	paired (TRUE), or any (NA) read should be returned.}

  \item{isProperPair}{A logical(1) indicating whether improperly paired
	(FALSE), properly paired (TRUE), or any (NA) read should be
	returned. A properly paired read is defined by the alignment
	algorithm and might, e.g., represent reads aligning to identical
	reference sequences and with a specified distance.}

  \item{isUnmappedQuery}{A logical(1) indicating whether unmapped
  (FALSE), mapped (TRUE), or any (NA) read should be returned.}

  \item{hasUnmappedMate}{A logical(1) indicating whether reads with
  mapped (FALSE), unmapped (TRUE), or any (NA) mate should be returned.}

  \item{isMinusStrand}{A logical(1) indicating whether reads aligned to
	the plus (FALSE), minus (TRUE), or any (NA) strand should be
	returned.}

  \item{isMateMinusStrand}{A logical(1) indicating whether mate reads
	aligned to the plus (FALSE), minus (TRUE), or any (NA) strand should
	be returned.}

  \item{isFirstMateRead}{A logical(1) indicating whether the first mate
	read should be returned (TRUE) or not (FALSE), or whether mate read
	number should be ignored (NA).}

  \item{isSecondMateRead}{A logical(1) indicating whether the second mate
	read should be returned (TRUE) or not (FALSE), or whether mate read
	number should be ignored (NA).}

  \item{isPrimaryRead}{A logical(1) indicating whether reads that are
	not primary (FALSE), are primary (TRUE) or whose primary status does
	not matter (NA) should be returned. A non-primary read might result
	when portions of a read aligns to multiple locations, e.g., when
	spanning splice junctions).}

  \item{isValidVendorRead}{A logical(1) indicating whether invalid
	(FALSE), valid (TRUE), or any (NA) read should be returned. A
	'valid' read is one flagged by the vendor as passing quality control
	criteria.}

  \item{isDuplicate}{A logical(1) indicating that un-duplicated (FALSE),
	duplicated (TRUE), or any (NA) reads should be
	returned. 'Duplicated' reads may represent PCR or optical
	duplicates.}

  \item{object}{An instance of class \code{ScanBamParam}.}

}

\section{Objects from the Class}{

  Objects are created by calls of the form \code{ScanBamParam()}.

}
\section{Slots}{
  \describe{

    \item{\code{flag}}{Object of class \code{integer} encoding flags
      to be kept when they have their '0' (\code{keep0}) or '1'
      (\code{keep1}) bit set.}

    \item{\code{simpleCigar}}{Object of class \code{logical}
     indicating, when TRUE, that only 'simple' cigars (empty or 'M') are
     returned.}

    \item{\code{which}}{Object of class \code{RangesList} indicating
      which reference sequence and coordinate reads must overlap.}

    \item{\code{what}}{Object of class \code{character} indicating
    which fields are to be returned.}

  }
}

\section{Functions and methods}{

  See 'Usage' for details on invocation.

  Constructor:
  \describe{

	\item{ScanBamParam:}{Returns a \code{ScanBamParam} object.}

  }

  Accessors:
  \describe{

	\item{bamWhat}{Returns a \code{character} vector of fields to be
	extracted.}

	\item{bamWhich}{Returns a \code{RangesList} of bounds on reads to
	  be extracted. A length 0 \code{RangesList} represents all reads.}

	\item{bamFlag}{Returns an \code{integer(2)} representation of reads
	  flagged to be kept or excluded.}

	\item{bamSimpleCigar}{Returns a \code{logical(1)} vector indicating
	  whether reads without indels or clipping be kept.}

  }

  Methods:
  \describe{

	\item{show}{Compactly display the object.}

  }

}
\author{Martin Morgan}

\seealso{
  \code{\link{scanBam}}
}
\examples{
# defaults
p0 <- ScanBamParam()

# subset of reads based on genomic coordinates
which <- RangesList(seq1=IRanges(1000, 2000),
                    seq2=IRanges(c(100, 1000), c(1000, 2000)))
p1 <- ScanBamParam(which=which)

# subset of reads based on 'flag' value
p2 <- ScanBamParam(flag=scanBamFlag(isMinusStrand=FALSE))

# subset of fields
p3 <- ScanBamParam(what=c("rname", "strand", "pos", "width"))
                
# use
fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
res <- scanBam(fl, param=p2)[[1]]
lapply(res, head)

}
\keyword{classes}
