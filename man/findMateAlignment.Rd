\name{findMateAlignment}

\alias{findMateAlignment}
\alias{findMateAlignment2}
\alias{makeGappedAlignmentPairs}

\alias{getDumpedAlignments}
\alias{countDumpedAlignments}
\alias{flushDumpedAlignments}

\title{Pairing the elements of a GappedAlignments object}

\description{
Utilities for pairing the elements of a \link[GenomicRanges]{GappedAlignments}
object.
}

\usage{
findMateAlignment(x, verbose=FALSE)
makeGappedAlignmentPairs(x, use.names=FALSE, use.mcols=FALSE)

## Related low-level utilities:
getDumpedAlignments()
countDumpedAlignments()
flushDumpedAlignments()
}

\arguments{
  \item{x}{
    A named \link[GenomicRanges]{GappedAlignments} object with metadata
    columns \code{flag}, \code{mrnm}, and \code{mpos}. Typically obtained
    by loading aligned paired-end reads from a BAM file with:
    \preformatted{
    param <- ScanBamParam(what=c("flag", "mrnm", "mpos"))
    x <- readBamGappedAlignments(..., use.names=TRUE, param=param)
    }
  }
  \item{verbose}{
    If \code{TRUE}, then \code{findMateAlignment} will print some
    details about what is currently going on. Mostly useful for debugging.
  }
  \item{use.names}{
    Whether the names on the input object should be propagated to the
    returned object or not.
  }
  \item{use.mcols}{
    Names of the metadata columns to propagate to the returned
    \link[GenomicRanges]{GappedAlignmentPairs} object.
  }
}

\details{
\subsection{Pairing algorithm used by findMateAlignment}{
\code{findMateAlignment} is the power horse used by higher-level
functions like \code{makeGappedAlignmentPairs} and
\code{\link{readBamGappedAlignmentPairs}} for pairing the records
loaded from a BAM file containing aligned paired-end reads.

It implements the following pairing algorithm:
\itemize{
  \item First only records with flag bit 0x1 set to 1, flag bit 0x4 set to 0,
        and flag bit 0x8 set to 0 are candidates for pairing (see the SAM
        Spec for a description of flag bits and fields).
        \code{findMateAlignment} will ignore any other record. That is,
        records that correspond to single-end reads, and records that
        correspond to paired-end reads where one or both ends are unmapped,
        are discarded.

  \item Then the algorithm looks at the following fields and flag bits:
        \itemize{
          \item (A) QNAME
          \item (B) RNAME, RNEXT
          \item (C) POS, PNEXT
          \item (D) Flag bits Ox10 and 0x20
          \item (E) Flag bits 0x40 and 0x80
          \item (F) Flag bit 0x2
          \item (G) Flag bit 0x100
        }
        2 records rec(i) and rec(j) are considered mates iff all the following
        conditions are satisfied:
        \itemize{
          \item (A) They have the same QNAME
          \item (B) RNEXT(i) == RNAME(j) and RNEXT(j) == RNAME(i)
          \item (C) PNEXT(i) == POS(j) and PNEXT(j) == POS(i)
          \item (D) Flag bit 0x20 of rec(i) == Flag bit 0x10 of rec(j) and
                    Flag bit 0x20 of rec(j) == Flag bit 0x10 of rec(i)
          \item (E) rec(i) corresponds to the first segment in the template and
                    rec(j) corresponds to the last segment in the template, OR,
                    rec(j) corresponds to the first segment in the template and
                    rec(i) corresponds to the last segment in the template
          \item (F) rec(i) and rec(j) have same flag bit 0x2
          \item (G) rec(i) and rec(j) have same flag bit 0x100
        }
}
}

\subsection{Ambiguous pairing}{
The above algorithm will find almost all pairs unambiguously, even when
the same pair of reads maps to several places in the genome. Note
that, when a given pair maps to a single place in the genome, looking
at (A) is enough to pair the 2 corresponding records. The additional
conditions (B), (C), (D), (E), (F), and (G), are only here to help in
the situation where more than 2 records share the same QNAME. And that
works most of the times. Unfortunately there are still situations where
this is not enough to solve the pairing problem unambiguously.

For example, here are 4 records (loaded in a GappedAlignments object)
that cannot be paired with the above algorithm:

Showing the 4 records as a GappedAlignments object of length 4:
\preformatted{
GappedAlignments with 4 alignments and 2 metadata columns:
                    seqnames strand       cigar    qwidth     start    end
                       <Rle>  <Rle> <character> <integer> <integer> <integer>
  SRR031714.2658602    chr2R      +  21M384N16M        37   6983850 6984270
  SRR031714.2658602    chr2R      +  21M384N16M        37   6983850 6984270
  SRR031714.2658602    chr2R      -  13M372N24M        37   6983858 6984266
  SRR031714.2658602    chr2R      -  13M378N24M        37   6983858 6984272
                        width      ngap |     mrnm      mpos
                    <integer> <integer> | <factor> <integer>
  SRR031714.2658602       421         1 |    chr2R   6983858
  SRR031714.2658602       421         1 |    chr2R   6983858
  SRR031714.2658602       409         1 |    chr2R   6983850
  SRR031714.2658602       415         1 |    chr2R   6983850
}
Note that the BAM fields show up in the following columns:
\itemize{
  \item QNAME: the names of the GappedAlignments object (unnamed col)
  \item RNAME: the seqnames col
  \item POS: the start col
  \item RNEXT: the mrnm col
  \item PNEXT: the mpos col
}
As you can see, the aligner has aligned the same pair to the same
location twice! The only difference between the 2 aligned pairs is in
the CIGAR i.e. one end of the pair is aligned twice to the same location
with exactly the same CIGAR while the other end of the pair is aligned
twice to the same location but with slightly different CIGARs.

Now showing the corresponding flag bits:
\preformatted{
     isPaired isProperPair isUnmappedQuery hasUnmappedMate isMinusStrand
[1,]        1            1               0               0             0
[2,]        1            1               0               0             0
[3,]        1            1               0               0             1
[4,]        1            1               0               0             1
     isMateMinusStrand isFirstMateRead isSecondMateRead isNotPrimaryRead
[1,]                 1               0                1                0
[2,]                 1               0                1                0
[3,]                 0               1                0                0
[4,]                 0               1                0                0
     isNotPassingQualityControls isDuplicate
[1,]                           0           0
[2,]                           0           0
[3,]                           0           0
[4,]                           0           0
}

As you can see, rec(1) and rec(2) are second mates, rec(3) and rec(4)
are both first mates. But looking at (A), (B), (C), (D), (E), (F), and (G),
the pairs could be rec(1) <-> rec(3) and rec(2) <-> rec(4), or they could
be rec(1) <-> rec(4) and rec(2) <-> rec(3). There is no way to
disambiguate!

So \code{findMateAlignment} is just ignoring (with a warning) those alignments
with ambiguous pairing, and dumping them in a place from which they can be
retrieved later (i.e. after \code{findMateAlignment} has returned) for
further examination (see "Dumped alignments" subsection below for the details).
In other words, alignments that cannot be paired unambiguously are not paired
at all. Concretely, this means that \code{\link{readGappedAlignmentPairs}}
is guaranteed to return a \link[GenomicRanges]{GappedAlignmentPairs} object
where every pair was formed in an non-ambiguous way. Note that, in practice,
this approach doesn't seem to leave aside a lot of records because ambiguous
pairing events seem pretty rare.
}

\subsection{Dumped alignments}{
Alignments with ambiguous pairing are dumped in a place ("the dump
environment") from which they can be retrieved with
\code{getDumpedAlignments()} after \code{findMateAlignment} has returned.

Two additional utilities are provided for manipulation of the dumped
alignments: \code{countDumpedAlignments} for counting them (a fast equivalent
to \code{length(getDumpedAlignments())}), and \code{flushDumpedAlignments} to
flush "the dump environment". Note that "the dump environment" is
automatically flushed at the beginning of a call to \code{findMateAlignment}.
}
}

\value{
  For \code{findMateAlignment}: An integer vector of the same length as
  \code{x}, containing only positive or NA values, where the i-th element
  is interpreted as follow:
  \itemize{
    \item An NA value means that no mate or more than 1 mate was found for
          \code{x[i]}.
    \item A non-NA value j gives the index in \code{x} of \code{x[i]}'s mate.
  }

  For \code{makeGappedAlignmentPairs}: A
  \link[GenomicRanges]{GappedAlignmentPairs} object where the pairs
  are formed internally by calling \code{findMateAlignment} on \code{x}.

  For \code{getDumpedAlignments}: \code{NULL} or a
  \link[GenomicRanges]{GappedAlignments} object containing the dumped
  alignments. See "Dumped alignments" subsection in the "Details" section
  above for the details.

  For \code{countDumpedAlignments}: The number of dumped alignments.

  Nothing for \code{flushDumpedAlignments}.
}

\author{H. Pages}

\seealso{
  \link[GenomicRanges]{GappedAlignments-class},
  \link[GenomicRanges]{GappedAlignmentPairs-class},
  \code{\link{readBamGappedAlignments}},
  \code{\link{readBamGappedAlignmentPairs}}
}

\examples{
bamfile <- system.file("extdata", "ex1.bam", package="Rsamtools")
param <- ScanBamParam(what=c("flag", "mrnm", "mpos"))
x <- readBamGappedAlignments(bamfile, use.names=TRUE, param=param)
mate <- findMateAlignment(x)
head(mate)
table(is.na(mate))
galp0 <- makeGappedAlignmentPairs(x)
galp <- makeGappedAlignmentPairs(x, use.name=TRUE, use.mcols="flag")
galp
colnames(mcols(galp))
colnames(mcols(first(galp)))
colnames(mcols(last(galp)))
}

\keyword{manip}
