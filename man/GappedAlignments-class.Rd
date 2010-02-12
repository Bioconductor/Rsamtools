\name{GappedAlignments-class}
\Rdversion{1.1}

\alias{class:GappedAlignments}
\alias{GappedAlignments-class}
\alias{GappedAlignments}
\alias{class:Alignments0}
\alias{Alignments0-class}
\alias{Alignments0}
\alias{class:Alignments1}
\alias{Alignments1-class}
\alias{Alignments1}

\alias{length,Alignments0-method}
\alias{rname}
\alias{rname,Alignments0-method}
\alias{rname,Alignments1-method}
\alias{strand,Alignments0-method}
\alias{strand,Alignments1-method}
\alias{cigar}
\alias{cigar,Alignments0-method}
\alias{cigar,Alignments1-method}
\alias{ranges,Alignments0-method}
\alias{ranges,Alignments1-method}
\alias{qwidth}
\alias{qwidth,GappedAlignments-method}
\alias{as.data.frame,GappedAlignments-method}
\alias{show,GappedAlignments-method}
\alias{readBAMasAlignments0}
\alias{readBAMasAlignments0,character-method}
\alias{readBAMasAlignments1}
\alias{readBAMasAlignments1,character-method}
\alias{[,Alignments0,ANY,ANY-method}
\alias{[,Alignments1,ANY,ANY-method}
\alias{coverage,GappedAlignments-method}


\title{GappedAlignments objects}

\description{
  The GappedAlignments class is a simple virtual container which
  purpose is to store a set of alignments that will hold just enough
  information for supporting the operations described below.

  There are currently 2 concrete implementations of the container:
  Alignments0 and Alignments1. Both implementations are equivalent
  from a user point of view i.e. they hold the same data (but organized
  differently) and allow the same operations. This is a temporary
  situation until the developpers make a choice.

  WARNING! This is work-in-progress. Expect frequent changes in
  functionalities. The name of the containers are provisory and
  reflecting the fact that we are at a very early stage of their
  design/implementation.
}

\details{
  A GappedAlignments object is a vector-like object where each element
  describes an alignment i.e. how a given sequence (called "query"
  or "read", typically short) aligns to a reference sequence (typically
  long).

  Most of the time, a GappedAlignments object will be created by loading
  records from a BAM (or SAM) file and each element in the resulting
  object will correspond to a record. BAM/SAM records generally contain
  a lot of information but only part of that information is loaded
  in the GappedAlignments object. In particular, we discard the query
  sequences (SEQ field), the query ids (QNAME field), the query qualities
  (QUAL), the mapping qualities (MAPQ) and any other information that
  is not needed in order to support the operations or methods described
  below.

  This means that multi-reads (i.e. reads with multiple hits in the
  reference) won't receive any special treatment i.e. the various SAM/BAM
  records corresponding to a multi-read will show up in the GappedAlignments
  object as if they were coming from different/unrelated queries.
  Also paired-end reads will be treated as single-end reads and the
  pairing information will be lost.

  Each element of a GappedAlignments object consists of:
  \itemize{
    \item The name of the reference sequence. (This is the RNAME field
          in a SAM/BAM record.)
    \item The strand in the reference sequence to which the query is
          aligned. (This information is stored in the FLAG field in a
          SAM/BAM record.)
    \item The CIGAR string in the "Extended CIGAR format" (see the SAM
          Format Specifications for the details).
    \item The 1-based leftmost position/coordinate of the clipped query
          relative to the reference sequence. We will refer to it as
          the "start" of the query. (This is the POS field in a SAM/BAM
          record.)
    \item The 1-based rightmost position/coordinate of the clipped query
          relative to the reference sequence. We will refer to it as
          the "end" of the query. (This is NOT explicitly stored in a
          SAM/BAM record but can be inferred from the POS and CIGAR fields.)
          Note that all positions/coordinates are always relative to
          the first base at the 5' end of the plus strand of the reference
          sequence, even when the query is aligned to the minus strand.
    \item The genomic intervals between the "start" and "end" of the query
          that are "covered" by the alignment. Saying that the full
          [start,end] interval is covered is the same as saying that the
          alignment has no gap (no N in the CIGAR). It is then considered
          a simple alignment. Note that a simple alignment can have
          mismatches or deletions (in the reference). In other words, a
          deletion, encoded with a D, is NOT considered a gap.
  }

  Strictly speaking, the last 2 items can be inferred from the CIGAR
  and the "start" so, in theory, they don't need to be explicitely
  stored in the GappedAlignments object. However, both implementations
  (Alignments0 and Alignments1) store it for performance reasons but
  at the cost of some redundancy (which impact on memory footprint
  seems to be acceptable).
  
  The rest of this man page will focus on describing how to:
  \itemize{
    \item Access the information stored in a GappedAlignments object
          in a way that is independent from how the data are actually
          stored internally.
    \item How to create and manipulate a GappedAlignments object.
  }
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{readBAMasAlignments0(file, index, ..., which)},
      \code{readBAMasAlignments1(file, index, ..., which)}:
      Loads a BAM file into an Alignments0 (or Alignments1) object.
      See \code{?\link{scanBam}} for a description of the arguments.
      Unlike SAM/BAM records, we don't support unaligned queries so
      we discard those records.
    }
  }
}

\section{Accessor methods}{
  In the code snippets below, \code{x} is a GappedAlignments object.

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
      \code{qwidth(x)}:
      Returns an integer vector of length \code{length(x)}
      containing the length of the query *after* hard clipping
      (i.e. the length of the query sequence that is stored in
      the corresponding SAM/BAM record).
    }
    \item{}{
      \code{ranges(x)}:
      Returns a \link[IRanges]{CompressedNormalIRangesList} object of
      length \code{length(x)} where each element represents regions in
      the reference to which a query is aligned. See Details section
      above for more information.
    }
    \item{}{
      \code{start(x)}, \code{end(x)}:
      Returns an integer vector of length \code{length(x)}
      containing the "start" and "end" (respectively) of the query
      for each alignment. See Details section above for the exact
      definitions of the "start" and "end" of a query.
      Note that \code{start(x)} and \code{end(x)} are equivalent
      to \code{min(ranges(x))} and \code{max(ranges(x))},
      respectively.
    }
    \item{}{
      \code{width(x)}:
      Defined as \code{end(x) - start(x) + 1L}.
      Note that this is generally different from \code{qwidth(x)}
      except for alignments with a trivial CIGAR string (i.e. a
      string of the form \code{"<n>M"} where <n> is a number).
    }
    \item{}{
      \code{ngap(x)}:
      Returns an integer vector of length \code{length(x)}
      containing the number of gaps for each alignment.
      Equivalent to \code{elementLengths(ranges(x)) - 1L}.
    }
  }
}

\section{Subsetting and related operations}{
  In the code snippets below, \code{x} is a GappedAlignments object.

  \describe{
    \item{}{
      \code{x[i]}:
      Returns a new GappedAlignments object (of the same type as \code{x},
      i.e. Alignments0 or Alignments1) made of the selected alignments.
      \code{i} can be a numeric or logical vector.
    }
  }
}

\section{Other methods}{
  In the code snippets below, \code{x} is a GappedAlignments object.

  \describe{
    \item{}{
      \code{coverage(x)}:
      Returns a named \link[IRanges]{RleList} object with one element
      (integer-Rle) per unique reference sequence. Each element represents
      \code{x}'s coverage of the corresponding reference sequence, that is,
      how many times each nucleotide position in the sequence is covered
      by the alignments in \code{x}.
      Note that the semantic of the \code{coverage} method for
      GappedAlignments objects is different from the semantic of the method
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
al0 <- readBAMasAlignments0(f1)
al0

## ---------------------------------------------------------------------
## A. BASIC MANIPULATION
## ---------------------------------------------------------------------
length(al0)
head(al0)
head(rname(al0))
head(strand(al0))
head(cigar(al0))
head(qwidth(al0))
table(qwidth(al0))

head(ranges(al0))
head(start(al0))
head(end(al0))
head(width(al0))
head(ngap(al0))

## ---------------------------------------------------------------------
## B. SUBSETTING
## ---------------------------------------------------------------------
al0[strand(al0) == "-"]
al0[grep("I", cigar(al0), fixed=TRUE)]
al0[grep("N", cigar(al0), fixed=TRUE)]  # no gaps

## A confirmation that all the queries map to the reference with no
## gaps:
stopifnot(all(elementLengths(ranges(al0)) == 1))
ranges(al0)[[1]]
ranges(al0)[[2]]
ii <- grep("D", cigar(al0), fixed=TRUE)
al0[ii]
## Ds are NOT gaps:
ranges(al0)[[ii[1]]]
ranges(al0)[[ii[2]]]

al0[qwidth(al0) != width(al0)]

## This MUST return an empty object:
al0[cigar(al0) == "35M" & qwidth(al0) != 35]
## but this doesn't have too:
al0[cigar(al0) != "35M" & qwidth(al0) == 35]



## ---------------------------------------------------------------------
## C. COVERAGE
## ---------------------------------------------------------------------
coverage(al0)
}

\keyword{methods}
\keyword{classes}
