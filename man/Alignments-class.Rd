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
\alias{qwidth}
\alias{qwidth,Alignments0-method}
\alias{gappedRanges}
\alias{gappedRanges,Alignments0-method}
\alias{start,Alignments0-method}
\alias{end,Alignments0-method}
\alias{as.data.frame,Alignments0-method}
\alias{show,Alignments0-method}
\alias{readBAMasAligments}
\alias{coerce,Alignments0,GappedRanges-method}
\alias{coerce,Alignments0,CompressedIRangesList-method}
\alias{coerce,Alignments0,IRangesList-method}
\alias{coerce,Alignments0,RangesList-method}
\alias{[,Alignments0,ANY,ANY-method}
\alias{[[,Alignments0,ANY,ANY-method}
\alias{elementLengths,Alignments0-method}
\alias{coverage,Alignments0-method}


\title{Alignments0 objects}

\description{
  The Alignments0 class is a simple container for storing a set of
  alignments that will hold just enough information for supporting
  the operations described below.

  WARNING! This is work-in-progess. Expect frequent changes in
  functionalities. The name of the container is provisory and
  reflecting the fact that we are at a very early stage of its
  design/implementation.
}

\details{
  An Alignments0 object is a vector-like object where each element
  describes an alignment i.e. how a given sequence (called "query"
  or "read", typically short) aligns to a reference sequence (typically
  long).

  Most of the times, an Alignments0 object will be created by loading
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
    \item The 1-based leftmost position/coordinate of the clipped query
          relative to the reference sequence. We will refer to it as
          the "start" of the query. (This is the POS field in a SAM/BAM
          record.)
    \item The 1-based rightmost position/coordinate of the clipped query
          relative to the reference sequence. We will refer to it as
          the "end" of the query. (This is NOT explicitly stored in a
          SAM/BAM record but can be inferred from the POS and CIGAR fields.)
    \item The 1-based positions/coordinates of the "gaps" in the
          reference (aka "skipped regions from the reference"). Note
          that a gap is encoded with an N in the CIGAR (a deletion,
          encoded with a D, is NOT considered a gap).
  }
  Strictly speaking, the last 2 items can be inferred from the CIGAR
  and the "start" so, in theory, they don't need to be explicitely
  stored in the Alignments0 object. In fact, the internals of the
  Alignments0 container could change and we don't want to commit to
  any particular internal representation.

  Hence the rest of this man page will focus on describing how to:
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
      \code{readBAMasAligments(file, index=file, which=RangesList())}:
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
      \code{qwidth(x)}:
      Returns an integer vector of length \code{length(x)}
      containing the length of the query *after* hard clipping
      (i.e. the length of the query sequence that is stored in
      the corresponding SAM/BAM record).
    }
    \item{}{
      \code{start(x)}, \code{end(x)}:
      Returns an integer vector of length \code{length(x)}
      containing the "start" and "end" (respectively) of the query
      for each alignment. See Details section above for the exact
      definitions of the "start" and "end" of a query.
    }
    \item{}{
      \code{width(x)}:
      Defined as \code{end(x) - start(x) + 1L}.
      Note that this is generally different from \code{qwidth(x)}
      except for alignments with a trivial CIGAR string (i.e. a
      string of the form \code{"<n>M"} where <n> is a number).
    }
    \item{}{
      \code{gappedRanges(x)}:
      Returns a \link[IRanges]{GappedRanges} object of length
      \code{length(x)} where each element is a
      \link[IRanges]{NormalIRanges} object representing a "gapped range".
      More on this below.
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
    \item{}{
      \code{x[[i]]}:
      This is a convenience for \code{gappedRanges(x)[[i]]}.
      Both return the \link[IRanges]{NormalIRanges} object
      containing the regions in the reference (1 range per region)
      to which the i-th query is aligned.
      \code{i} must be a single integer.
    }
    \item{}{
      \code{elemenType(x)}:
      Returns the type of \code{x[[i]]}.
      Note that the semantic of the \code{[[} method for Alignments0
      objects is different from the semantic of the method for
      \link[IRanges]{Ranges} objects (the latter returns an integer
      vector).
    }
    \item{}{
      \code{elementLengths(x)}:
      This is a convenience for \code{elementLengths(gappedRanges(x))}.
      Returns an integer vector of length \code{length(x)}
      containing the number of ranges in x[[i]].
      Semantically equivalent to
      \preformatted{sapply(seq_len(length(x)), function(i) length(x[[i]]))}
      but much faster.
      Note that the semantic of the \code{elementLengths} method for
      Alignments0 objects is different from the semantic of the method
      for \link[IRanges]{Ranges} objects (the latter returns the
      \code{width} of the \link[IRanges]{Ranges} object).
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
  H. Pages
}

\seealso{
  \code{\link{scanBam}},
  \link[IRanges]{NormalIRanges-class},
  \link[IRanges]{GappedRanges-class},
  \code{\link[IRanges]{coverage}},
  \link[IRanges]{RleList-class}
}

\examples{
f1 <- system.file("extdata", "ex1.bam", package="Rsamtools")
al1 <- readBAMasAligments(f1)
al1

## ---------------------------------------------------------------------
## A. BASIC MANIPULATION
## ---------------------------------------------------------------------
length(al1)
head(al1)
head(rname(al1))
head(strand(al1))
head(cigar(al1))
head(qwidth(al1))
head(start(al1))
head(end(al1))
head(qwidth(al1))
table(qwidth(al1))
gappedRanges(al1)

## ---------------------------------------------------------------------
## B. SUBSETTING
## ---------------------------------------------------------------------
al1[strand(al1) == "-"]
al1[grep("I", cigar(al1), fixed=TRUE)]
al1[grep("N", cigar(al1), fixed=TRUE)]  # no gaps

## A confirmation that all the queries map to the reference with no
## gaps:
stopifnot(all(elementLengths(al1) == 1))
al1[[1]]
al1[[2]]
ii <- grep("D", cigar(al1), fixed=TRUE)
al1[ii]
## Ds are NOT gaps:
al1[[ii[1]]]
al1[[ii[2]]]

al1[qwidth(al1) != width(al1)]

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
