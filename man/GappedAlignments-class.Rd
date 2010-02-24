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

\alias{length,GappedAlignments-method}
\alias{rname}
\alias{rname,Alignments0-method}
\alias{rname,Alignments1-method}
\alias{rname<-}
\alias{rname<-,Alignments0-method}
\alias{rname<-,Alignments1-method}
\alias{strand,Alignments0-method}
\alias{strand,Alignments1-method}
\alias{cigar}
\alias{cigar,GappedAlignments-method}
\alias{granges,Alignments0-method}
\alias{granges,Alignments1-method}
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
\alias{coerce,Alignments1,Alignments0-method}
\alias{coerce,Alignments0,Alignments1-method}
\alias{[,Alignments0,ANY,ANY-method}
\alias{[,Alignments1,ANY,ANY-method}
\alias{shift,Alignments0-method}
\alias{shift,Alignments1-method}
\alias{updateCigarAndStart,Alignments0-method}
\alias{updateCigarAndStart,Alignments1-method}
\alias{qnarrow,GappedAlignments-method}
\alias{coverage,GappedAlignments-method}
\alias{findOverlaps,GappedAlignments,ANY-method}
\alias{findOverlaps,ANY,GappedAlignments-method}
\alias{findOverlaps,GappedAlignments,GappedAlignments-method}


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
      \code{rname(x) <- value}:
      Replace the name of the reference sequence for each alignment.
      \code{value} must be a character factor/vector, or a 'character' Rle,
      or a 'factor' Rle, with the same length as \code{x}.
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
      \code{granges(x)}, \code{ranges(x)}:
      Returns a \link[BSgenome]{GRangesList} or
      \link[IRanges]{CompressedNormalIRangesList} object of
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

  \describe{
    \item{}{
      \code{shift(x, shift)}:
      \code{x} is a GappedAlignments object.
      Returns a new GappedAlignments object (of the same type as \code{x},
      i.e. Alignments0 or Alignments1) where all the alignments have been
      shifted along the reference by the number of positions specified in
      the \code{shift} argument (vector of integers).
    }
    \item{}{
      \code{qnarrow(x, start=NA, end=NA, width=NA)}:
      \code{x} is a GappedAlignments object.
      Returns a new GappedAlignments object (of the same type as \code{x},
      i.e. Alignments0 or Alignments1) where all the query sequences have
      been narrowed according to the \code{start}, \code{end} and
      \code{width} arguments.
    }
    \item{}{
      \code{coverage(x)}:
      \code{x} is a GappedAlignments object.
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
    \item{}{
      \code{findOverlaps(query, subject, ...)}:
      \code{query} or \code{subject} or both are GappedAlignments objects.
      Equivalent to \code{findOverlaps(granges(query), subject, ...)}
      when \code{query} is a GappedAlignments object, or to
      \code{findOverlaps(query, granges(subject), ...)} when
      \code{subject} is a GappedAlignments object, or to
      \code{findOverlaps(granges(query), granges(subject), ...)}
      when both are GappedAlignments objects.
      See \link[BSgenome]{findOverlaps,GRangesList,GRanges-method}
      in the BSgenome package for more information (in particular for
      a description of the extra arguments and the returned object).
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
  \link[BSgenome]{GRangesList-class},
  \link[IRanges]{NormalIRanges-class},
  \link[IRanges]{CompressedNormalIRangesList-class},
  \code{\link[IRanges]{coverage}},
  \link[IRanges]{RleList-class},
  \link[BSgenome]{findOverlaps,GRangesList,GRanges-method}
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
levels(rname(al0))

## Rename the reference sequences:
rname(al0) <- sub("seq", "chr", rname(al0))
levels(rname(al0))

head(strand(al0))
head(cigar(al0))
head(qwidth(al0))
table(qwidth(al0))

granges(al0)  # a GRangesList object
ranges(al0)   # a CompressedNormalIRangesList object
stopifnot(identical(elementLengths(granges(al0)), elementLengths(ranges(al0))))

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
stopifnot(all(ngap(al0) == 0))

## Different ways to subset:
al0[6]             # a GappedAlignments object of length 1
granges(al0)[[6]]  # a GRanges object of length 1
ranges(al0)[[6]]   # a NormalIRanges object of length 1

## Ds are NOT gaps:
ii <- grep("D", cigar(al0), fixed=TRUE)
al0[ii]
ngap(al0[ii])
granges(al0[ii])

## qwidth() vs width():
al0[qwidth(al0) != width(al0)]

## This MUST return an empty object:
al0[cigar(al0) == "35M" & qwidth(al0) != 35]
## but this doesn't have too:
al0[cigar(al0) != "35M" & qwidth(al0) == 35]

## ---------------------------------------------------------------------
## C. qnarrow()
## ---------------------------------------------------------------------
## This trims 3 nucleotides on the left and 5 nucleotides on the right
## of each alignment:
qnarrow(al0, start=4, end=-6)
## Note that the 'start' and 'end' arguments specify what part of each
## query sequence should be kept (negative values being relative to the
## right end of the query sequence), not what part should be trimmed.

## Trimming on the left doesn't change the "end" of the queries.
qnarrow(al0, start=21)
stopifnot(identical(end(qnarrow(al0, start=21)), end(al0)))

## ---------------------------------------------------------------------
## D. coverage()
## ---------------------------------------------------------------------
coverage(al0)

## ---------------------------------------------------------------------
## E. findOverlaps()
## ---------------------------------------------------------------------
findOverlaps(al0, granges(al0)[[1]])
}

\keyword{methods}
\keyword{classes}
