\name{Alignments-class}
\Rdversion{1.1}

\alias{class:Alignments0}
\alias{Alignments0-class}
\alias{Alignments0}


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
  Alignments0 container could change and we don't want to comit to
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
      See \code{?\link{scabBam}} for a description of the arguments.
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
    }
    \item{}{
      \code{rname(x)}:
    }
    \item{}{
      \code{strand(x)}:
    }
    \item{}{
      \code{cigar(x)}:
    }
    \item{}{
      \code{qwidth(x)}:
    }
    \item{}{
      \code{start(x)}:
    }
    \item{}{
      \code{end(x)}:
    }
    \item{}{
      \code{width(x)}:
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
  \link[IRanges]{IRanges-class},
  \link[IRanges]{GappedRanges-class},
  \code{\link[IRanges]{coverage}},
  \link[IRanges]{RleList-class}
}

\examples{
}

\keyword{methods}
\keyword{classes}
