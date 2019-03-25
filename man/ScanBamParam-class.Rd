\name{ScanBamParam}
\Rdversion{1.1}
\docType{class}
\alias{ScanBamParam-class}
\alias{ScanBamParam}
\alias{ScanBamParam,missing-method}
\alias{ScanBamParam,ANY-method}
\alias{ScanBamParam,IntegerRangesList-method}
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
\alias{bamTagFilter}
\alias{bamTagFilter<-}
\alias{bamWhat<-}
\alias{bamWhat}
\alias{bamWhich<-}
\alias{bamWhich<-,ScanBamParam,GRanges-method}
\alias{bamWhich<-,ScanBamParam,IntegerRangesList-method}
\alias{bamWhich<-,ScanBamParam,ANY-method}
\alias{bamWhich}
\alias{bamMapqFilter}
\alias{bamMapqFilter<-}
% methods
\alias{show,ScanBamParam-method}
% flag utils
\alias{FLAG_BITNAMES}
\alias{bamFlagAsBitMatrix}
\alias{bamFlagAND}
\alias{bamFlagTest}
% flag names
\alias{isPaired}
\alias{isProperPair}
\alias{isUnmappedQuery}
\alias{hasUnmappedMate}
\alias{isMinusStrand}
\alias{isMateMinusStrand}
\alias{isFirstMateRead}
\alias{isSecondMateRead}
\alias{isNotPrimaryRead} % deprecated
\alias{isSecondaryAlignment}
\alias{isNotPassingQualityControls}
\alias{isDuplicate}

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
    reverseComplement = FALSE, tag = character(0), tagFilter = list(),
    what = character(0), which, mapqFilter=NA_integer_)

# Constructor helpers
scanBamFlag(isPaired = NA, isProperPair = NA, isUnmappedQuery = NA, 
    hasUnmappedMate = NA, isMinusStrand = NA, isMateMinusStrand = NA,
    isFirstMateRead = NA, isSecondMateRead = NA, isNotPrimaryRead = NA,
    isSecondaryAlignment = NA, isNotPassingQualityControls = NA,
    isDuplicate = NA, isSupplementaryAlignment = NA)

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
bamTagFilter(object)
bamTagFilter(object) <- value
bamWhat(object)
bamWhat(object) <- value
bamWhich(object)
bamWhich(object) <- value
bamMapqFilter(object)
bamMapqFilter(object) <- value

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

  \item{reverseComplement}{A logical(1) vectors. BAM files store reads
    mapping to the minus strand as though they are on the plus
    strand. Rsamtools obeys this convention by default
    (\code{reverseComplement=FALSE}), but when this value is set to TRUE
    returns the sequence and quality scores of reads mapped to the minus
    strand in the reverse complement (sequence) and reverse (quality) of
    the read as stored in the BAM file. This might be useful if wishing
    to recover read and quality scores as represented in fastq files,
    but is NOT appropriate for variant calling or other alignment-based
    operations.}

  \item{tag}{A character vector naming tags to be extracted. A tag is an
    optional field, with arbitrary information, stored with each
    record. Tags are identified by two-letter codes, so all elements of
    \code{tag} must have exactly 2 characters.}

  \item{tagFilter}{A named \code{list} of atomic vectors. The name of
    each \code{list} element is the tag name (two-letter code), and the
    corresponding atomic vector is the set of acceptable values for the
    tag. Only reads with specified tags are included. \code{NULL}s,
    \code{NA}s, and empty strings are not allowed in the atomic
    vectors.}

  \item{what}{A character vector naming the fields to return
    \code{scanBamWhat()} returns a vector of available fields.
    Fields are described on the \code{\link{scanBam}} help page.}

  \item{mapqFilter}{A non-negative integer(1) specifying the minimum
    mapping quality to include. BAM records with mapping qualities less
    than \code{mapqFilter} are discarded.}

  \item{which}{A \code{\linkS4class{GRanges}},
    \code{\linkS4class{IntegerRangesList}}, or any object that can be
    coerced to a \code{IntegerRangesList}, or missing object, from which
    a \code{IRangesList} instance will be constructed. Names of the
    \code{IRangesList} correspond to reference sequences, and ranges to
    the regions on that reference sequence for which matches are
    desired. Because data types are coerced to \code{IRangesList},
    \code{which} does \emph{not} include strand information (use the
    \code{flag} argument instead). Only records with a read overlapping
    the specified ranges are returned. All ranges must have ends less
    than or equal to 536870912. When one record overlaps two ranges in
    \code{which}, the record is returned \emph{twice}.}

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

  \item{isNotPrimaryRead}{Deprecated; use \code{isSecondaryAlignment}.}

  \item{isSecondaryAlignment}{A logical(1) indicating whether alignments
    that are secondary (TRUE), are not (FALSE) or whose secondary status
    does not matter (NA) should be returned. A non-primary alignment
    (\dQuote{secondary alignment} in the SAM specification) might result
    when a read aligns to multiple locations. One alignment is
    designated as primary and has this flag set to FALSE; the remainder,
    for which this flag is TRUE, are designated by the aligner as
    secondary.}

  \item{isNotPassingQualityControls}{A logical(1) indicating whether
    reads passing quality controls (FALSE), reads not passing quality
    controls (TRUE), or any (NA) read should be returned.}

  \item{isDuplicate}{A logical(1) indicating that un-duplicated (FALSE),
    duplicated (TRUE), or any (NA) reads should be
    returned. 'Duplicated' reads may represent PCR or optical
    duplicates.}

  \item{isSupplementaryAlignment}{A logical(1) indicating that
    non-supplementary (FALSE), supplementary (TRUE), or any (NA) reads
    should be returned. The SAM specification indicates that
    'supplementary' reads are part of a chimeric alignment.}
  
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

    \item{\code{tagFilter}}{Object of class \code{list} (named)
      indicating tags to filter by, and the set of acceptable values for
      each tag.}

    \item{\code{what}}{Object of class \code{character} indicating
      what fields are to be returned.}

    \item{\code{which}}{Object of class \code{IntegerRangesList} indicating
      which reference sequence and coordinate reads must overlap.}

    \item{\code{mapqFilter}}{Object of class \code{integer} indicating
      the minimum mapping quality required for input, or NA to indicate
      no filtering.}

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

    \item{bamTagFilter, bamTagFilter<-}{Returns or sets a named
      \code{list} of tags to filter by, and the set of their acceptable
      values.}
  
    \item{bamWhat, bamWhat<-}{Returns or sets a \code{character} vector
      of fields to be extracted.}

    \item{bamWhich, bamWhich<-}{Returns or sets a \code{IntegerRangesList} of
      bounds on reads to be extracted. A length 0 \code{IntegerRangesList}
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
which <- IRangesList(seq1=IRanges(1000, 2000),
                     seq2=IRanges(c(100, 1000), c(1000, 2000)))
p1 <- ScanBamParam(what=scanBamWhat(), which=which)

## subset of reads based on 'flag' value
p2 <- ScanBamParam(what=scanBamWhat(),
                   flag=scanBamFlag(isMinusStrand=FALSE))

## subset of fields
p3 <- ScanBamParam(what=c("rname", "strand", "pos", "qwidth"))
                
## use
fl <- system.file("extdata", "ex1.bam", package="Rsamtools",
                  mustWork=TRUE)
res <- scanBam(fl, param=p2)[[1]]
lapply(res, head)

## tags; NM: edit distance; H1: 1-difference hits
p4 <- ScanBamParam(tag=c("NM", "H1"), what="flag")
bam4 <- scanBam(fl, param=p4)
str(bam4[[1]][["tag"]])

## tagFilter
p5 <- ScanBamParam(tag=c("NM", "H1"), tagFilter=list(NM=c(2, 3, 4)))
bam5 <- scanBam(fl, param=p5)
table(bam5[[1]][["tag"]][["NM"]])

## flag utils
flag <- scanBamFlag(isUnmappedQuery=FALSE, isMinusStrand=TRUE)

p6 <- ScanBamParam(what="flag")
bam6 <- scanBam(fl, param=p6)
flag6 <- bam6[[1]][["flag"]]
head(bamFlagAsBitMatrix(flag6[1:9]))
colSums(bamFlagAsBitMatrix(flag6))
flag
bamFlagAsBitMatrix(flag)
}
\keyword{classes}
