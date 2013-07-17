\name{sequenceLayer}

\alias{sequenceLayer}


\title{Lay read sequences alongside the reference space, using their CIGARs}

\description{
  \code{sequenceLayer} can lay strings that belong to a given space (e.g.
  the \code{"query"} space) alongside another space (e.g. the
  \code{"reference"} space) by removing/injecting substrings from/into them,
  using the supplied CIGARs.

  Its primary use case is to lay the read sequences stored in a BAM file
  (which are considered to belong to the \code{"query"} space) alongside
  the \code{"reference"} space. It can also be used to remove the parts
  of the read sequences that correspond to soft-clipping. More generally
  it can lay strings that belong to any supported space alongside any other
  supported space. See the Details section below for the list of supported
  spaces.
}

\usage{
sequenceLayer(x, cigar, from="query", to="reference",
              D.letter="-", N.letter="-",
              I.letter="-", S.letter="+", H.letter="+")
}

\arguments{
  \item{x}{
    An \link{XStringSet} object containing strings that belong to a given
    space.
  }
  \item{cigar}{
    A character vector or factor of the same length as \code{x} containing
    the extended CIGAR strings (one per element in \code{x}).
  }
  \item{from, to}{
    A single string specifying one of the 8 supported spaces listed in the
    Details section below. \code{from} must be the current space (i.e. the
    space the strings in \code{x} belong to) and \code{to} is the space
    alonside which to lay the strings in \code{x}.
  }
  \item{D.letter, N.letter, I.letter, S.letter, H.letter}{
    A single letter used as a filler for injections. More on this in
    the Details section below.
  }
}

\details{
  The 8 supported spaces are: \code{"reference"},
  \code{"reference-N-regions-removed"}, \code{"query"},
  \code{"query-before-hard-clipping"}, \code{"query-after-soft-clipping"},
  \code{"pairwise"}, \code{"pairwise-N-regions-removed"},
  and \code{"pairwise-dense"}.

  Each space can be characterized by the extended CIGAR
  operations that are \emph{visible} in it. A CIGAR operation is
  said to be \emph{visible} in a given space if it "runs along it",
  that is, if it's associated with a block of contiguous positions in
  that space (the size of the block being the length of the operation).
  For example, the M/=/X operations are \emph{visible} in all spaces,
  the D/N operations are \emph{visible} in the \code{"reference"} space
  but not in the \code{"query"} space, the S operation is \emph{visible}
  in the \code{"query"} space but not in the \code{"reference"}
  or in the \code{"query-after-soft-clipping"} space, etc...

  Here are the extended CIGAR operations that are \emph{visible} in each
  space:
  \enumerate{
    \item reference: M, D, N, =, X

    \item reference-N-regions-removed: M, D, =, X

    \item query: M, I, S, =, X

    \item query-before-hard-clipping: M, I, S, H, =, X

    \item query-after-soft-clipping: M, I, =, X

    \item pairwise: M, I, D, N, =, X

    \item pairwise-N-regions-removed: M, I, D, =, X

    \item pairwise-dense: M, =, X
  }

  \code{sequenceLayer} lays a string that belongs to one space alongside
  another by (1) removing the substrings associated with operations that
  are not \emph{visible} anymore in the new space, and (2) injecting
  substrings associated with operations that become \emph{visible} in the
  new space. Each injected substring has the length of the operation
  associated with it, and its content is controlled via the corresponding
  \code{*.letter} argument.

  For example, when going from the \code{"query"} space to the
  \code{"reference"} space (the default), the I- and S-substrings (i.e.
  the substrings associated with I/S operations) are removed, and
  substrings associated with D/N operations are injected. More precisely,
  the D-substrings are filled with the letter specified in \code{D.letter},
  and the N-substrings with the letter specified in \code{N.letter}.
  The other \code{*.letter} arguments are ignored in that case.
}

\value{
  An \link{XStringSet} object of the same class and length as \code{x}.
}

\author{H. Pages}

\seealso{
  \itemize{
    \item The \code{\link{stackStringsFromBam}} function
          for stacking the read sequences (or their quality strings)
          stored in a BAM file on a region of interest.

    \item The \code{\link{readGAlignmentsFromBam}} function for loading read
          sequences from a BAM file (via a \link[GenomicRanges]{GAlignments}
          object).

    \item The \code{\link[Biostrings]{extractAt}} and
          \code{\link[Biostrings]{replaceAt}} functions in the \pkg{Biostrings}
          package for extracting/replacing arbitrary substrings from/in a
          string or set of strings.

    \item \link[GenomicRanges]{cigar-utils} in the \pkg{GenomicRanges}
          package for the CIGAR utility functions used internally by
          \code{sequenceLayer}.
  }
}

\examples{
## ---------------------------------------------------------------------
## A. FROM "query" TO "reference" SPACE
## ---------------------------------------------------------------------

## Load read sequences from a BAM file (they will be returned in a
## GAlignments object):
bamfile <- system.file("extdata", "ex1.bam", package="Rsamtools")
param <- ScanBamParam(what="seq")
gal <- readGAlignmentsFromBam(bamfile, param=param)
qseq <- mcols(gal)$seq  # the read sequences (aka query sequences)

## Lay the query sequences alongside the reference space. This will
## remove the substrings associated with insertions to the reference
## (I operations) and soft clipping (S operations), and will inject new
## substrings (filled with "-") where deletions from the reference (D
## operations) and skipped regions from the reference (N operations)
## occurred during the alignment process:
qseq_on_ref <- sequenceLayer(qseq, cigar(gal))

## A typical use case for doing the above is to compute 1 consensus
## sequence per chromosome. The code below shows how this can be done
## in 2 extra steps.

## Step 1: Compute one consensus matrix per chromosome.
qseq_on_ref_by_chrom <- splitAsList(qseq_on_ref, seqnames(gal))
pos_by_chrom <- splitAsList(start(gal), seqnames(gal))

cm_by_chrom <- lapply(names(pos_by_chrom),
    function(seqname)
        consensusMatrix(qseq_on_ref_by_chrom[[seqname]],
                        as.prob=TRUE,
                        shift=pos_by_chrom[[seqname]]-1,
                        width=seqlengths(gal)[[seqname]]))
names(cm_by_chrom) <- names(pos_by_chrom)

## 'cm_by_chrom' is a list of consensus matrices. Each matrix has 17
## rows (1 per letter in the DNA alphabet) and 1 column per chromosome
## position.

## Step 2: Compute the consensus string from each consensus matrix.
## We'll put "+" in the strings wherever there is no coverage for that
## position, and "N" where there is coverage but no consensus.
cs_by_chrom <- lapply(cm_by_chrom,
    function(cm) {
        ## Because consensusString() doesn't like consensus matrices
        ## with columns that contain only zeroes (and you will have
        ## columns like that for chromosome positions that don't
        ## receive any coverage), we need to "fix" 'cm' first.
        idx <- colSums(cm) == 0
        cm["+", idx] <- 1
        DNAString(consensusString(cm, ambiguityMap="N"))
    })

## consensusString() provides some flexibility to let you extract
## the consensus in different ways. See '?consensusString' in the
## Biostrings package for the details.

## Finally, note that the read quality strings can also be used as
## input for sequenceLayer():
param <- ScanBamParam(what="qual")
gal <- readGAlignmentsFromBam(bamfile, param=param)
qual <- mcols(gal)$qual  # the read quality strings
qual_on_ref <- sequenceLayer(qual, cigar(gal))
## Note that since the "-" letter is a valid quality code, there is
## no way to distinguish it from the "-" letters inserted by
## sequenceLayer().

## ---------------------------------------------------------------------
## B. FROM "query" TO "query-after-soft-clipping" SPACE
## ---------------------------------------------------------------------

## Going from "query" to "query-after-soft-clipping" simply removes
## the substrings associated with soft clipping (S operations):
qseq <- DNAStringSet(c("AAAGTTCGAA", "TTACGATTAN", "GGATAATTTT"))
cigar <- c("3H10M", "2S7M1S2H", "2M1I1M3D2M4S")
clipped_qseq <- sequenceLayer(qseq, cigar,
                              from="query", to="query-after-soft-clipping")

sequenceLayer(clipped_qseq, cigar,
              from="query-after-soft-clipping", to="query")

sequenceLayer(clipped_qseq, cigar,
              from="query-after-soft-clipping", to="query",
              S.letter="-")

## ---------------------------------------------------------------------
## C. BRING QUERY AND REFERENCE SEQUENCES TO THE "pairwise" or
##    "pairwise-dense" SPACE
## ---------------------------------------------------------------------

## Load read sequences from a BAM file:
library(RNAseqData.HNRNPC.bam.chr14)
bamfile <- RNAseqData.HNRNPC.bam.chr14_BAMFILES[1]
param <- ScanBamParam(what="seq")
gal <- readGAlignmentsFromBam(bamfile, param=param)
qseq <- mcols(gal)$seq  # the read sequences (aka query sequences)

## Load the corresponding reference sequences from the appropriate
## BSgenome package (the reads in RNAseqData.HNRNPC.bam.chr14 were
## aligned to hg19):
library(BSgenome.Hsapiens.UCSC.hg19)
rseq <- getSeq(Hsapiens, as(gal, "GRanges"))  # the reference sequences

## Bring 'qseq' and 'rseq' to the "pairwise" space.
## For 'qseq', this will remove the substrings associated with soft
## clipping (S operations) and inject substrings (filled with "-")
## associated with deletions from the reference (D operations) and
## skipped regions from the reference (N operations). For 'rseq', this
## will inject substrings (filled with "-") associated with insertions
## to the reference (I operations).
qseq2 <- sequenceLayer(qseq, cigar(gal),
                       from="query", to="pairwise")
rseq2 <- sequenceLayer(rseq, cigar(gal),
                       from="reference", to="pairwise")

## Sanity check: 'qseq2' and 'rseq2' should have the same shape.
stopifnot(identical(elementLengths(qseq2), elementLengths(rseq2)))

## A closer look at reads with insertions and deletions:
cigar_op_table <- cigarOpTable(cigar(gal))
head(cigar_op_table)

I_idx <- which(cigar_op_table[ , "I"] >= 2)  # at least 2 insertions
qseq2[I_idx]
rseq2[I_idx]

D_idx <- which(cigar_op_table[ , "D"] >= 2)  # at least 2 deletions
qseq2[D_idx]
rseq2[D_idx]

## A closer look at reads with skipped regions:
N_idx <- which(cigar_op_table[ , "N"] != 0)
qseq2[N_idx]
rseq2[N_idx]

## A variant of the "pairwise" space is the "pairwise-dense" space.
## In that space, all indels and skipped regions are removed from 'qseq'
## and 'rseq'.
qseq3 <- sequenceLayer(qseq, cigar(gal),
                       from="query", to="pairwise-dense")
rseq3 <- sequenceLayer(rseq, cigar(gal),
                       from="reference", to="pairwise-dense")

## Sanity check: 'qseq3' and 'rseq3' should have the same shape.
stopifnot(identical(elementLengths(qseq3), elementLengths(rseq3)))

## Insertions were removed:
qseq3[I_idx]
rseq3[I_idx]

## Deletions were removed:
qseq3[D_idx]
rseq3[D_idx]

## Skipped regions were removed:
qseq3[N_idx]
rseq3[N_idx]

## ---------------------------------------------------------------------
## D. SANITY CHECKS
## ---------------------------------------------------------------------
SPACES <- c("reference",
            "reference-N-regions-removed",
            "query",
            "query-before-hard-clipping",
            "query-after-soft-clipping",
            "pairwise",
            "pairwise-N-regions-removed",
            "pairwise-dense")

cigarWidth <- list(
    function(cigar) cigarWidthAlongReferenceSpace(cigar),
    function(cigar) cigarWidthAlongReferenceSpace(cigar,
                                                  N.regions.removed=TRUE),
    function(cigar) cigarWidthAlongQuerySpace(cigar),
    function(cigar) cigarWidthAlongQuerySpace(cigar,
                                              before.hard.clipping=TRUE),
    function(cigar) cigarWidthAlongQuerySpace(cigar,
                                              after.soft.clipping=TRUE),
    function(cigar) cigarWidthAlongPairwiseSpace(cigar),
    function(cigar) cigarWidthAlongPairwiseSpace(cigar,
                                                 N.regions.removed=TRUE),
    function(cigar) cigarWidthAlongPairwiseSpace(cigar, dense=TRUE)
)

cigar <- c("3H2S4M1D2M2I1M5N3M6H", "5M1I3M2D4M2S")

seq <- list(
    BStringSet(c(A="AAAA-BBC+++++DDD", B="AAAAABBB--CCCC")),
    BStringSet(c(A="AAAA-BBCDDD", B="AAAAABBB--CCCC")),
    BStringSet(c(A="..AAAABBiiCDDD", B="AAAAAiBBBCCCC..")),
    BStringSet(c(A=".....AAAABBiiCDDD......", B="AAAAAiBBBCCCC..")),
    BStringSet(c(A="AAAABBiiCDDD", B="AAAAAiBBBCCCC")),
    BStringSet(c(A="AAAA-BBiiC+++++DDD", B="AAAAAiBBB--CCCC")),
    BStringSet(c(A="AAAA-BBiiCDDD", B="AAAAAiBBB--CCCC")),
    BStringSet(c(A="AAAABBCDDD", B="AAAAABBBCCCC"))
)

stopifnot(all(sapply(1:8,
    function(i) identical(width(seq[[i]]), cigarWidth[[i]](cigar))
)))

sequenceLayer2 <- function(x, cigar, from, to)
    sequenceLayer(x, cigar, from=from, to=to, D.letter="-", N.letter= "+",
                  I.letter="i", S.letter=".", H.letter=".")

identical_XStringSet <- function(target, current)
{
    ok1 <- identical(class(target), class(current))
    ok2 <- identical(names(target), names(current))
    ok3 <- all(target == current)
    ok1 && ok2 && ok3
}

res <- sapply(1:8, function(i) {
           sapply(1:8, function(j) {
               target <- seq[[j]]
               current <- sequenceLayer2(seq[[i]], cigar,
                                         from=SPACES[i], to=SPACES[j])
               identical_XStringSet(target, current)
           })
       })
stopifnot(all(res))
}

\keyword{methods}
\keyword{manip}
