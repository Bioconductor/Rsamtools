\name{ScanBamParam}
\Rdversion{1.1}
\docType{class}
\alias{ScanBamParam-class}
\alias{ScanBamParam}
\alias{ScanBamParam,missing-method}
\alias{ScanBamParam,RangesList-method}
\alias{ScanBamParam,RangedData-method}
\alias{ScanBamParam,GRanges-method}
% helpers
\alias{scanBamWhat}
\alias{scanBamFlag}
% accessors
\alias{bamFlag<-}
\alias{bamFlag}
\alias{bamReverseComplement<-}
\alias{bamReverseComplement}
\alias{bamSimpleCigar<-}
\alias{bamSimpleCigar}
\alias{bamTag<-}
\alias{bamTag}
\alias{bamWhat<-}
\alias{bamWhat}
\alias{bamWhich<-}
\alias{bamWhich<-,ScanBamParam,GRanges-method}
\alias{bamWhich<-,ScanBamParam,RangedData-method}
\alias{bamWhich<-,ScanBamParam,RangesList-method}
\alias{bamWhich<-,ScanBamParam,ANY-method}
\alias{bamWhich}
% methods
\alias{show,ScanBamParam-method}
% flag utils
\alias{FLAG_BITNAMES}
\alias{bamFlagAsBitMatrix}
\alias{bamFlagAND}
\alias{bamFlagTest}

\title{Parameters for scanning BAM files}
\description{

  Use \code{ScanBamParam()} to create a parameter object influencing
  what fields and which records are imported from a (binary) BAM file.
  Use of \code{which} requires that a BAM index file
  (\code{<filename>.bai}) exists.

}

\usage{

# Constructor
ScanBamParam(flag = scanBamFlag(), simpleCigar = FALSE,
    reverseComplement = FALSE, tag = character(0),
    what = character(0), which)

# Constructor helpers
scanBamFlag(isPaired = NA, isProperPair = NA, isUnmappedQuery = NA, 
    hasUnmappedMate = NA, isMinusStrand = NA, isMateMinusStrand = NA,
    isFirstMateRead = NA, isSecondMateRead = NA, isNotPrimaryRead = NA,
    isNotPassingQualityControls = NA, isDuplicate = NA,
    isValidVendorRead = NA)

scanBamWhat()

# Accessors
bamFlag(object, asInteger=FALSE)
bamFlag(object) <- value
bamReverseComplement(object)
bamReverseComplement(object) <- value
bamSimpleCigar(object)
bamSimpleCigar(object) <- value
bamTag(object)
bamTag(object) <- value
bamWhat(object)
bamWhat(object) <- value
bamWhich(object)
bamWhich(object) <- value

\S4method{show}{ScanBamParam}(object)

# Flag utils
bamFlagAsBitMatrix(flag, bitnames=FLAG_BITNAMES)
bamFlagAND(flag1, flag2)
bamFlagTest(flag, value)
}

\arguments{

  \item{flag}{For \code{ScanBamParam}, an integer(2) vector used to filter
    reads based on their 'flag' entry. This is most easily created with the
    \code{scanBamFlag()} helper function.

    For \code{bamFlagAsBitMatrix}, \code{bamFlagTest} an integer vector
    where each element represents a 'flag' entry.  }

  \item{simpleCigar}{A logical(1) vector which, when TRUE, returns only
    those reads for which the cigar (run-length encoded representation
    of the alignment) is missing or contains only matches / mismatches
    (\code{'M'}).}

  \item{reverseComplement}{A logical(1) vector which, when TRUE, returns
    the sequence and quality scores of reads mapped to the minus strand
    in the reverse complement (sequence) and reverse (quality) of the
    read as stored in the BAM file.}

  \item{tag}{A character vector naming tags to be extracted. A tag is an
    optional field, with arbitrary information, stored with each
    record. Tags are identified by two-letter codes, so all elements of
    \code{tag} must have exactly 2 characters.}

  \item{what}{A character vector naming the fields to
    return. \code{scanBamWhat()} returns a vector of available
    fields. Fields are described on the \code{\link{scanBam}} help
    page.}

  \item{which}{A \code{\linkS4class{GRanges}},
    \code{\linkS4class{RangesList}}, \code{\linkS4class{RangedData}}, or
    missing object, from which a \code{IRangesList} instance will be
    constructed. Names of the \code{IRangesList} correspond to reference
    sequences, and ranges to the regions on that reference sequence for
    which matches are desired. Because data types are coerced to
    \code{IRangesList}, \code{which} does \emph{not} include strand
    information (use the \code{flag} argument instead). Only records
    with a read overlapping the specified ranges are returned. All
    ranges must have ends less than or equal to 536870912.}

  \item{isPaired}{A logical(1) indicating whether unpaired (FALSE),
    paired (TRUE), or any (NA) read should be returned.}

  \item{isProperPair}{A logical(1) indicating whether improperly paired
    (FALSE), properly paired (TRUE), or any (NA) read should be
    returned. A properly paired read is defined by the alignment
    algorithm and might, e.g., represent reads aligning to identical
    reference sequences and with a specified distance.}

  \item{isUnmappedQuery}{A logical(1) indicating whether unmapped
    (TRUE), mapped (FALSE), or any (NA) read should be returned.}

  \item{hasUnmappedMate}{A logical(1) indicating whether reads with
    mapped (FALSE), unmapped (TRUE), or any (NA) mate should be
    returned.}

  \item{isMinusStrand}{A logical(1) indicating whether reads aligned to
    the plus (FALSE), minus (TRUE), or any (NA) strand should be
    returned.}

  \item{isMateMinusStrand}{A logical(1) indicating whether mate reads
    aligned to the plus (FALSE), minus (TRUE), or any (NA) strand should
    be returned.}

  \item{isFirstMateRead}{A logical(1) indicating whether the first mate
    read should be returned (TRUE) or not (FALSE), or whether mate read
    number should be ignored (NA).}

  \item{isSecondMateRead}{A logical(1) indicating whether the second
    mate read should be returned (TRUE) or not (FALSE), or whether mate
    read number should be ignored (NA).}

  \item{isNotPrimaryRead}{A logical(1) indicating whether reads that are
    primary (FALSE), are not primary (TRUE) or whose primary status does
    not matter (NA) should be returned. A non-primary read might result
    when portions of a read aligns to multiple locations, e.g., when
    spanning splice junctions).}

  \item{isNotPassingQualityControls}{A logical(1) indicating whether
    reads passing quality controls (FALSE), reads not passing quality
    controls (TRUE), or any (NA) read should be returned.}

  \item{isValidVendorRead}{Deprecated; use
    \code{isNotPassingQualityControls}.}

  \item{isDuplicate}{A logical(1) indicating that un-duplicated (FALSE),
    duplicated (TRUE), or any (NA) reads should be
    returned. 'Duplicated' reads may represent PCR or optical
    duplicates.}
  
  \item{object}{An instance of class \code{ScanBamParam}.}

  \item{value}{An instance of the corresponding slot, to be assigned to
    \code{object} or, for \code{bamFlagTest}, a \code{character(1)} name
    of the flag to test, e.g., \dQuote{isUnmappedQuery}, from the
    arguments to \code{scanBamFlag}.}

  \item{asInteger}{logical(1) indicating whether \sQuote{flag} should be
    returned as an encoded integer vector (\code{TRUE}) or
    human-readable form (\code{FALSE}).}

  \item{bitnames}{Names of the flag bits to extract. Will be the colnames
    of the returned matrix.}

  \item{flag1, flag2}{Integer vectors containing \sQuote{flag} entries.}
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

    \item{\code{reverseComplement}}{Object of class \code{logical}
      indicating, when TRUE, that reads on the minus strand are to be
      reverse complemented (sequence) and reversed (quality).}

    \item{\code{tag}}{Object of class \code{character} indicating what
      tags are to be returned.}

    \item{\code{what}}{Object of class \code{character} indicating
      what fields are to be returned.}

    \item{\code{which}}{Object of class \code{RangesList} indicating
      which reference sequence and coordinate reads must overlap.}

  }
}

\section{Functions and methods}{

  See 'Usage' for details on invocation.

  Constructor:
  \describe{

    \item{ScanBamParam:}{Returns a \code{ScanBamParam} object. The
      \code{which} argument to the constructor can be one of several
      different types, as documented above.}

  }

  Accessors:
  \describe{

    \item{bamTag, bamTag<-}{Returns or sets a \code{character} vector of
      tags to be extracted.}

    \item{bamWhat, bamWhat<-}{Returns or sets a \code{character} vector
      of fields to be extracted.}

    \item{bamWhich, bamWhich<-}{Returns or sets a \code{RangesList} of
      bounds on reads to be extracted. A length 0 \code{RangesList}
      represents all reads.}

    \item{bamFlag, bamFlag<-}{Returns or sets an \code{integer(2)}
      representation of reads flagged to be kept or excluded.}

    \item{bamSimpleCigar, bamSimpleCigar<-}{Returns or sets a
      \code{logical(1)} vector indicating whether reads without indels
      or clipping be kept.}

    \item{bamReverseComplement, bamReverseComplement<-}{Returns or sets
      a \code{logical(1)} vector indicating whether reads on the minus
      strand will be returned with sequence reverse complemented and
      quality reversed.}

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
## defaults
p0 <- ScanBamParam()

## subset of reads based on genomic coordinates
which <- RangesList(seq1=IRanges(1000, 2000),
                    seq2=IRanges(c(100, 1000), c(1000, 2000)))
p1 <- ScanBamParam(which=which)

## subset of reads based on 'flag' value
p2 <- ScanBamParam(flag=scanBamFlag(isMinusStrand=FALSE))

## subset of fields
p3 <- ScanBamParam(what=c("rname", "strand", "pos", "qwidth"))
                
## use
fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
res <- scanBam(fl, param=p2)[[1]]
lapply(res, head)

## tags; NM: edit distance; H1: 1-difference hits
p4 <- ScanBamParam(tag=c("NM", "H1"), what="flag")
bam4 <- scanBam(fl, param=p4)
str(bam4[[1]][["tag"]])

## flag utils
flag <- scanBamFlag(isUnmappedQuery=FALSE, isMinusStrand=TRUE)
flag
bamFlagAsBitMatrix(flag)
flag4 <- bam4[[1]][["flag"]]
bamFlagAsBitMatrix(flag4[1:9], bitnames=c("isUnmappedQuery", "isMinusStrand"))
}
\keyword{classes}
