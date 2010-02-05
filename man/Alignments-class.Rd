\name{Alignments-class}
\Rdversion{1.1}

\alias{class:Alignments0}
\alias{Alignments0-class}
\alias{Alignments0}

\alias{length,Alignments0-method}
\alias{rname}
\alias{rname,Alignments0-method}
\alias{strand,Alignments0-method}
\alias{cigar}
\alias{cigar,Alignments0-method}
\alias{ranges,Alignments0-method}
\alias{qwidth}
\alias{qwidth,Alignments0-method}
\alias{isSimple}
\alias{isSimple,Alignments0-method}
\alias{as.data.frame,Alignments0-method}
\alias{show,Alignments0-method}
\alias{readBAMasAlignments}
\alias{readBAMasAlignments,character-method}
\alias{[,Alignments0,ANY,ANY-method}
\alias{coverage,Alignments0-method}


\title{Alignments0 objects}

\description{
  The Alignments0 class is a simple container for storing a set of
  alignments that will hold just enough information for supporting
  the operations described below.

  WARNING! This is work-in-progress. Expect frequent changes in
  functionalities. The name of the container is provisory and
  reflecting the fact that we are at a very early stage of its
  design/implementation.
}

\details{
  An Alignments0 object is a vector-like object where each element
  describes an alignment i.e. how a given sequence (called "query"
  or "read", typically short) aligns to a reference sequence (typically
  long).

  Most of the time, an Alignments0 object will be created by loading
  records from a BAM (or SAM) file and each element in the resulting
  object will correspond to a record. BAM/SAM records generally contain
  a lot of information but only part of that information is loaded
  in the Alignments0 object. In particular, we discard the query
  sequences (SEQ field), the query ids (QNAME field), the query qualities
  (QUAL), the mapping qualities (MAPQ) and any other information that
  is not needed in order to support the operations or methods described
  below.

  This means that multi-reads (i.e. reads with multiple hits in the
  reference) won't receive any special treatment i.e. the various SAM/BAM
  records corresponding to a multi-read will show up in the Alignments0
  object as if they were coming from different/unrelated queries.
  Also paired-end reads will be treated as single-end reads and the
  pairing information will be lost.

  Each element of an Alignments0 object consists of:
  \itemize{
    \item The name of the reference sequence. (This is the RNAME field
          in a SAM/BAM record.)
    \item The strand in the reference sequence to which the query is
          aligned. (This information is stored in the FLAG field in a
          SAM/BAM record.)
    \item The CIGAR string in the "Extended CIGAR format" (see the SAM
          Format Specifications for the details).
    \item The 1-based intervals, i.e. [start, end] ranges, of the
          aligned clipped query relative to the reference sequence. A
          simple alignment is represented by a single interval, while a
          gapped alignment, which is defined to have non-aligned skipped
          regions (or gaps) on the reference, contains an interval for
          each region that is aligned. In CIGAR notation a non-aligned
          gap is coded as an 'N', which is not to be confused with the
          (aligned) deletion code 'D'.
  }

  The rest of this man page will focus on describing how to:
  \itemize{
    \item Access the information stored in an Alignments0 object in
          a way that is independent from how the data are actually
          stored internally.
    \item How to create and manipulate Alignments0 objects.
  }
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{readBAMasAlignments(file, index, ..., which)}:
      Loads a BAM file into an Alignments0 object.
      See \code{?\link{scanBam}} for a description of the arguments.
      Unlike SAM/BAM records, we don't support unaligned queries so
      we discard those records.
    }
  }
}

\section{Accessor methods}{
  In the code snippets below, \code{x} is an Alignments0 object.

  \describe{
    \item{}{
      \code{length(x)}:
      Returns the number of alignments in \code{x}.
    }
    \item{}{
      \code{rname(x)}:
      Returns a character factor of length \code{length(x)}
      containing the name of the reference sequence for each alignment.
    }
    \item{}{
      \code{strand(x)}:
      Returns a character factor of length \code{length(x)}
      (with levels +, - and *) containing the strand in the reference
      sequence to which the query is aligned.
    }
    \item{}{
      \code{cigar(x)}:
      Returns a character vector of length \code{length(x)}
      containing the CIGAR string for each alignment.
    }
    \item{}{
      \code{ranges(x)}:
      Returns a \link[IRanges]{CompressedNormalIRangesList} object of
      length \code{length(x)} where each element represents regions in
      the reference to which a query is aligned. See Details section
      above for more information.
    }
    \item{}{
      \code{qwidth(x)}:
      Returns an integer vector of length \code{length(x)}
      containing the length of the query *after* hard clipping
      (i.e. the length of the query sequence that is stored in
      the corresponding SAM/BAM record).
    }
    \item{}{
      \code{isSimple(x)}:
      Returns a logical vector of length \code{length(x)} whose
      elements indicate whether or not an alignment is simple
      (i.e. contains only matches and mismatches and is represented
      by a CIGAR string with a single 'M' run).
    }
  }
}

\section{Subsetting and related operations}{
  In the code snippets below, \code{x} is an Alignments0 object.

  \describe{
    \item{}{
      \code{x[i]}:
      Returns a new Alignments0 object made of the selected alignments.
      \code{i} can be a numeric or logical vector.
    }
  }
}

\section{Other methods}{
  In the code snippets below, \code{x} is an Alignments0 object.

  \describe{
    \item{}{
      \code{coverage(x)}:
      Returns a named \link[IRanges]{RleList} object with one element
      (integer-Rle) per unique reference sequence. Each element represents
      \code{x}'s coverage of the corresponding reference sequence, that is,
      how many times each nucleotide position in the sequence is covered
      by the alignments in \code{x}.
      Note that the semantic of the \code{coverage} method for
      Alignments0 objects is different from the semantic of the method
      for \link[IRanges]{Ranges} objects (the latter returns a single
      integer-Rle object representing the coverage of all ranges
      relatively to a unique imaginary reference sequence).
    }
  }
}

\references{
  \url{http://samtools.sourceforge.net/}
}

\author{
  H. Pages and P. Aboyoun
}

\seealso{
  \code{\link{scanBam}},
  \link[IRanges]{NormalIRanges-class},
  \link[IRanges]{CompressedNormalIRangesList-class},
  \code{\link[IRanges]{coverage}},
  \link[IRanges]{RleList-class}
}

\examples{
f1 <- system.file("extdata", "ex1.bam", package="Rsamtools")
al1 <- readBAMasAlignments(f1)
al1

## ---------------------------------------------------------------------
## A. BASIC MANIPULATION
## ---------------------------------------------------------------------
length(al1)
head(al1)
head(rname(al1))
head(strand(al1))
head(cigar(al1))
head(ranges(al1))
head(min(ranges(al1)))
head(max(ranges(al1)))
head(qwidth(al1))
table(qwidth(al1))

## ---------------------------------------------------------------------
## B. SUBSETTING
## ---------------------------------------------------------------------
al1[strand(al1) == "-"]
al1[grep("I", cigar(al1), fixed=TRUE)]
al1[grep("N", cigar(al1), fixed=TRUE)]  # no gaps

## A confirmation that all the queries map to the reference with no
## gaps:
stopifnot(all(elementLengths(ranges(al1)) == 1))
ranges(al1)[[1]]
ranges(al1)[[2]]
ii <- grep("D", cigar(al1), fixed=TRUE)
al1[ii]
## Ds are NOT gaps:
ranges(al1)[[ii[1]]]
ranges(al1)[[ii[2]]]

al1[!isSimple(al1)]

## This MUST return an empty object:
al1[cigar(al1) == "35M" & qwidth(al1) != 35]
## but this doesn't have too:
al1[cigar(al1) != "35M" & qwidth(al1) == 35]



## ---------------------------------------------------------------------
## C. COVERAGE
## ---------------------------------------------------------------------
coverage(al1)
}

\keyword{methods}
\keyword{classes}
