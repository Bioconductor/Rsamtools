\name{pileup}
\Rdversion{1.1}
% PileupParam
\alias{PileupParam-class}
\alias{PileupParam}
\alias{show,PileupParam-method}
\alias{max_depth}
\alias{min_base_quality}
\alias{min_mapq}
\alias{min_nucleotide_depth}
\alias{min_minor_allele_depth}
\alias{distinguish_strands}
\alias{distinguish_nucleotides}
\alias{ignore_query_Ns}
\alias{include_deletions}
\alias{include_insertions}
\alias{cycle_bins}

% pileup
\alias{pileup}
\alias{pileup,character-method}
\alias{pileup,BamFile-method}

\title{

  Use filters and output formats to calculate pile-up statistics for a
  BAM file.

}
\description{

  \code{pileup} uses \code{PileupParam} and \code{ScanBamParam} objects
  to calculate pileup statistics for a BAM file. The result is a
  \code{data.frame} with columns summarizing counts of reads overlapping
  each genomic position, optionally differentiated on nucleotide,
  strand, and position within read.

}
\usage{

pileup(file, index=file, ..., scanBamParam=ScanBamParam(),
       pileupParam=PileupParam())

## PileupParam constructor
PileupParam(max_depth=250, min_base_quality=0, min_mapq=13,
    min_nucleotide_depth=1, min_minor_allele_depth=0,
    distinguish_strands=TRUE, distinguish_nucleotides=TRUE,
    ignore_query_Ns=TRUE, include_deletions=TRUE, include_insertions=FALSE,
    cycle_bins=numeric())
}

\arguments{

  \item{file}{
    character(1) or \code{\link{BamFile}}; BAM file path.
  }

  \item{index}{
    When \code{file} is character(1), an optional character(1) of BAM
    index file path; see \code{\link{scanBam}}.
  }

  \item{\dots}{Additional arguments, perhaps used by methods.}

  \item{scanBamParam}{An instance of \code{\link{ScanBamParam}}.}

  \item{pileupParam}{An instance of \code{\link{PileupParam}}.}

  %% args for PileupParam

  \item{max_depth}{integer(1); maximum number of overlapping alignments
    considered for each position in the pileup.}

  \item{min_base_quality}{integer(1); minimum \sQuote{QUAL} value for
    each nucleotide in an alignment.}

  \item{min_mapq}{integer(1); minimum \sQuote{MAPQ} value for an
    alignment to be included in pileup.}

  \item{min_nucleotide_depth}{integer(1); minimum number of occurrences
    of a nucleotide at a given position required for a particular nucleotide to
    appear in the result.}

  \item{min_minor_allele_depth}{integer(1); minimum number of
    occurrences of \emph{all} nucleotides other than the major allele at
    a given position required for a particular nucleotide to appear in
    the result.}

  \item{distinguish_strands}{logical(1); \code{TRUE} if result should
    differentiate between strands.}

  \item{distinguish_nucleotides}{logical(1); \code{TRUE} if result
    should differentiate between nucleotides.}

  \item{ignore_query_Ns}{logical(1); \code{TRUE} if non-determinate
    nucleotides in alignments should be excluded from the pileup.}

  \item{include_deletions}{logical(1); \code{TRUE} to include deletions
    in pileup.}

  \item{include_insertions}{logical(1); \code{TRUE} to include
    insertions in pileup.}

  \item{cycle_bins}{numeric; non-negative, unique positions within a
    read to delimit bins. At least two values are required so at least
    one range can be formed. numeric(0) (default) indicates no
    binning. Sorted order is not required. Floating-point values are
    coerced to \code{integer}.}

% better wording?: 'if _____ should disqualify the read/alignment from
% the pileup'
}

\details{

  \emph{NB}: Use of \code{pileup} assumes familiarity with
  \code{\link{ScanBamParam}}, and use of \code{cycle_bins} assumes
  familiarity with \code{\link{cut}}.

  \code{pileup} visits each position in the BAM file, subject to
  constraints implied by \code{which} and \code{flag} of
  \code{scanBamParam}. For a given position, all reads overlapping the
  position that are consistent with constraints dictated by \code{flag}
  of \code{scanBamParam} and \code{pileupParam} are used for counting.

  If no \code{which} argument is supplied to the \code{ScanBamParam},
  behavior reflects that of \code{scanBam}: the entire file is visited
  and counted.
  
  Use a \code{\link{yieldSize}} value when creating a \code{BamFile}
  instance to manage memory consumption when using pileup with large BAM
  files. There are some differences in how \code{pileup} behavior is
  affected when the \code{yieldSize} value is set on the BAM file. See
  more details below.
  
  Many of the parameters of the \code{pileupParam} interact. A simple
  illustration is \code{ignore_query_Ns} and
  \code{distinguish_nucleotides}, as mentioned in the
  \code{ignore_query_Ns} section.

  Parameters for \code{pileupParam} belong to two categories: parameters
  that affect only the filtering criteria (so-called
  \sQuote{behavior-only} policies), and parameters that affect
  filtering behavior and the schema of the results
  (\sQuote{behavior+structure} policies).

  %% behavior-only

  %% (other behavior-only arguments are self-evident or sam
  %% spcification-evident)

  %% behavior+structure

  \code{distinguish_nucleotides} and \code{distinguish_strands} when set
  to \code{TRUE} each add a column (\code{nucleotide} and \code{strand},
  respectively) to the resulting \code{data.frame}. If both are TRUE,
  each combination of \code{nucleotide} and \code{strand} at a given
  position is counted separately. Setting only one to \code{TRUE}
  behaves as expected; for example, if only \code{nucleotide} is set to
  \code{TRUE}, each nucleotide at a given position is counted
  separately, but the distinction of on which strand the nucleotide
  appears is ignored.

  \code{ignore_query_Ns} determines how ambiguous nucloetides are
  treated. By default, ambiguous nucleotides (typically \sQuote{N} in
  BAM files) in alignments are ignored. If \code{ignore_query_Ns} is
  \code{FALSE}, ambiguous nucleotides are included in counts; further,
  if \code{ignore_query_Ns} is \code{FALSE} and
  \code{distinguish_nucleotides} is \code{TRUE} the \sQuote{N}
  nucleotide value appears in the nucleotide column when a base at a
  given position is ambiguous.

  By default, deletions with respect to the reference genome to which
  the reads were aligned are included in the counts in a pileup. If
  \code{include_deletions} is \code{TRUE} and
  \code{distinguish_nucleotides} is \code{TRUE}, the nucleotide column
  uses a \sQuote{-} character to indicate a deletion at a position.

  The \sQuote{=} nucleotide code in the \code{SEQ} field (to mean
  \sQuote{identical to reference genome}) is supported by pileup, such
  that a match with the reference genome is counted separately in the
  results if \code{distinguish_nucleotides} is \code{TRUE}.

  CIGAR support: \code{pileup} handles the extended CIGAR format by only
  heeding operations that contribute to counts (\sQuote{M}, \sQuote{D},
  \sQuote{I}). If you are confused about how the different CIGAR
  operations interact, see the SAM versions of the BAM files used for
  \code{pileup} unit testing for a fairly comprehensive human-readable
  treatment.
  
  \itemize{
    \item{Deletions and clipping:

      The extended CIGAR allows a number of operations conceptually
      similar to a \sQuote{deletion} with respect to the reference
      genome, but offer more specialized meanings than a simple
      deletion. CIGAR \sQuote{N} operations (not to be confused with
      \sQuote{N} used for non-determinate bases in the \code{SEQ} field)
      indicate a large skipped region, \sQuote{S} a soft clip, and
      \sQuote{H} a hard clip. \code{pileup} ignores \sQuote{N},
      \sQuote{S}, and \sQuote{H} CIGAR operations; only counts of true
      deletions (\sQuote{D} in the CIGAR) can be included by setting
      \code{include_deletions} to \code{TRUE}.

    }

    \item{Insertions and padding:

      \code{pileup} can include counts of insertion operations by
      settting \code{include_insertions} to \code{TRUE}. Details about
      insertions are effectively truncated such that each insertion is
      reduced to a single \sQuote{+} nucleotide. Information about
      e.g. the nucleotide code and base quality within the insertion is
      not included. The position of an insertion is the same as the
      called base that precedes it.

      Because \sQuote{+} is used as a nucleotide code to represent
      insertions in \code{pileup}, counts of the \sQuote{+} nucleotide
      participate in voting for \code{min_nucleotide_depth} and
      \code{min_minor_allele_depth} functionality.

      \sQuote{P} CIGAR (padding) operations never affect counts.

    }
  }

  \code{cycle_bins} allows users to differentiate pileup counts based on
  arbitrary regions within a read. \code{pileup} relies on
  \code{\link{cut}} to derive bins, but limits input to non-negative
  numeric values (coerced to \code{integer}s), including \code{Inf}. If
  \code{cycle_bins} is set \code{pileup} automatically excludes bases
  outside the implicit outer range. There are two important factors to
  consider when specifying bins:
  
  \itemize{
      \item{Because \code{\link{cut}} produces ranges in the form
      (first,last], \sQuote{0} should be used to create a bin that
      includes the first position.}

      \item{In order to account for variable-length reads in BAM files,
      \code{Inf} should be used to create a bin that extends to the end
      of all reads.}
  } %% end subsection cycle_bins gotchas

  \code{pileup} obeys \code{\link{yieldSize}} on \code{BamFile} objects
  with some differences from how \code{scanBam} responds to
  \code{yieldSize}. Here are some points to keep in mind when using
  \code{pileup} in conjunction with \code{yieldSize}:

  \itemize{

      \item{\code{BamFile}s with a \code{yieldSize} value set, once
      opened and used with \code{pileup}, \emph{should not be used} with
      other functions that accept a \code{BamFile} instance as
      input. Create a new \code{BamFile} instance instead of trying to
      reuse.}

      \item{\code{pileup} only returns genomic positions for which all
      input has been processed. \code{pileup} will hold in reserve
      positions for which only partial data has been
      processed. Positions held in reserve will be returned upon
      subsequent calls to \code{pileup} when all the input for a given
      position has been processed.}
    
      \item{The correspondence between yieldSize and the number of rows
      in the \code{data.frame} returned from \code{pileup} is not
      1-to-1. \code{yieldSize} only limits the number of
      \emph{alignments processed} from the BAM file each time the file
      is read. Only a handful of alignments can lead to many distinct
      records in the result.}

      \item{Like \code{scanBam}, \code{pileup} uses an empty return
      object (a zero-row \code{data.frame}) to indicate end-of-input.}

      \item{Sometimes reading \code{yieldSize} records from the BAM file
      does not result in any completed positions. In order to avoid
      returning a zero-row \code{data.frame} \code{pileup} will
      repeatedly process \code{yieldSize} additional records until at
      least one position can be returned to the user.}

  }

}

\value{
  
  For \code{pileup} a \code{data.frame} with 1 row per unique
  combination of differentiating column values that satisfied filter
  criteria, with frequency (\code{count}) of unique combination. Columns
  \code{seqnames}, \code{pos}, and \code{count} always appear; optional
  \code{strand}, \code{nulceotide}, and \code{cycle_bin} columns appear
  as dictated by arguments to \code{PileupParam}.

  Column details:

  \itemize{
    
    \item{\code{seqnames}: factor. SAM \sQuote{RNAME} field of
    alignment.}

    \item{\code{pos}: integer(1). Genomic position of base. Derived by
    offset from SAM \sQuote{POS} field of alignment.}

    \item{\code{strand}: factor. \sQuote{strand} to which read aligns.}

    \item{\code{nucleotide}: factor. \sQuote{nucleotide} of base,
      extracted from SAM \sQuote{SEQ} field of alignment.}

    \item{\code{cycle_bin}: factor. Cycle bin in which base appears.}

    \item{\code{count}: integer(1). Frequency of combination of
    differentiating fields, as indicated by values of other columns.}
    
  }

  See \code{\link{scanBam}} for more detailed explanation of SAM fields.

  If a \code{which} argument is specified for the \code{scanBamParam}, a
  \code{which_label} column (\code{factor} in the form
  \sQuote{rname:start-end}) is automatically included. The
  \code{which_label} column is used to maintain grouping of results,
  such that two queries of the same genomic region can be
  differentiated.

  Order of rows in \code{data.frame} is first by order of
  \code{seqnames} column based on the BAM index file, then
  non-decreasing order on columns \code{pos}, then \code{nucleotide},
  then \code{strand}, then \code{cycle_bin}.

  \code{PileupParam} returns an instance of PileupParam class.

}

\seealso{
  \itemize{
    \item{\link{Rsamtools}}
    \item{\linkS4class{BamFile}}
    \item{\link{ScanBamParam}}
    \item{\link{scanBam}}
    \item{\link{cut}}
  }
}

\author{Nate Hayden \url{nhayden@fhcrc.org}}

%% This information must go somewhere--essential for an outside
%% developer to help build context quickly
%% \excruciatingdetails{

%%   A careful reading of the C source code reveals that there are three
%%   stages of filtering; (1) filtering based on \code{scanBamParam}
%%   arguments, (2) filtering based on arguments of \code{pileupParam}
%%   that apply on a per-read and per-alignment basis and lastly, (3) a
%%   final filtering stage applied on the intermediate count of a
%%   position after all overlapping alignments have been
%%   processed. \code{min_nucleotide_depth} and
%%   \{min_minor_allele_depth} arguments cannot be applied until the
%%   last stage.

%% }

\examples{

## The examples below apply equally to pileup queries for specific
## genomic ranges (specified by the 'which' parameter of 'ScanBamParam')
## and queries across entire files; the entire-file examples are
## included first to make them easy to find. The more detailed examples
## of pileup use queries with a 'which' argument.

library("RNAseqData.HNRNPC.bam.chr14")

fl <- RNAseqData.HNRNPC.bam.chr14_BAMFILES[1]

\donttest{
## no 'which' argument to ScanBamParam: process entire file at once
res <- pileup(fl)
head(res)
table(res$strand, res$nucleotide)
}

\donttest{
## no 'which' argument to ScanBamParam with 'yieldSize' set on BamFile
bf <- open(BamFile(fl, yieldSize=5e4))
repeat {
    res <- pileup(bf)
    message(nrow(res), " rows in result data.frame")
    if(nrow(res) == 0L)
        break
}
close(bf)
}

## pileup for the first half of chr14 with default Pileup parameters
## 'which' argument: process only specified genomic range(s)
sbp <- ScanBamParam(which=GRanges("chr14", IRanges(1, 53674770)))
res <- pileup(fl, scanBamParam=sbp)
head(res)
table(res$strand, res$nucleotide)

## specify pileup parameters: include ambiguious nucleotides
## (the 'N' level in the nucleotide column of the data.frame)
p_param <- PileupParam(ignore_query_Ns=FALSE)
res <- pileup(fl, scanBamParam=sbp, pileupParam=p_param)
head(res)
table(res$strand, res$nucleotide)

## Don't distinguish strand, require a minimum frequency of 7 for a
## nucleotide at a genomic position to be included in results.

p_param <- PileupParam(distinguish_strands=TRUE,
                       min_nucleotide_depth=7)
res <- pileup(fl, scanBamParam=sbp, pileupParam=p_param)
head(res)
table(res$nucleotide, res$strand)

## Any combination of the filtering criteria is possible: let's say we
## want a "coverage pileup" that only counts reads with mapping
## quality >= 13, and bases with quality >= 10 that appear at least 4
## times at each genomic position
p_param <- PileupParam(distinguish_nucleotides=FALSE,
                       distinguish_strands=FALSE,
                       min_mapq=13,
                       min_base_quality=10,
                       min_nucleotide_depth=4)
res <- pileup(fl, scanBamParam=sbp, pileupParam=p_param)
head(res)
res <- res[, c("pos", "count")] ## drop seqnames and which_label cols
plot(count ~ pos, res, pch=".")

## Well-supported polymorphic positions (minor allele present in at
## least 5 alignments) with high map quality
p_param <- PileupParam(min_minor_allele_depth=5,
                       min_mapq=40,
                       distinguish_strand=FALSE)
res <- pileup(fl, scanBamParam=sbp, pileupParam=p_param)
dim(res) ## reduced to our biologically interesting positions
head(xtabs(count ~ pos + nucleotide, res))

## cycle_bins

\donttest{
## basic use: single pivot; count bases that appear in first 10
## positions of a read separately from the rest
p_param <- PileupParam(cycle_bins=c(0, 10, Inf))
res <- pileup(fl, scanBamParam=sbp, pileupParam=p_param)

## basic use: simple range; only include bases in positions 6-10 within
## read
p_param <- PileupParam(cycle_bins=c(5, 10))
res <- pileup(fl, scanBamParam=sbp, pileupParam=p_param)
}

## typical use: beginning, middle, and end segments; because of the
## nature of sequencing technology, it is common for bases in the
## beginning and end segments of each read to be less reliable. pileup
## makes it easy to count each segment separately.

## Assume determined ahead of time that for the data 1-7 makes sense for
## beginning, 8-12 middle, >=13 end (actual positions should be tailored
## to data in actual BAM files).
p_param <- PileupParam(cycle_bins=c(0, 7, 12, Inf)) ## note non-symmetric
res <- pileup(fl, scanBamParam=sbp, pileupParam=p_param)
xt <- xtabs(count ~ nucleotide + cycle_bin, res)
print(xt)
t(t(xt) / colSums(xt)) ## cheap normalization for illustrative purposes

## 'implicit outer range': in contrast to c(0, 7, 12, Inf), suppose we
##  still want to have beginning, middle, and end segements, but know
##  that the first three positions and any bases beyond the 16th
##  position are irrelevant. Hence, the implicit outer range is (3,16];
##  all bases outside of that are dropped.
p_param <- PileupParam(cycle_bins=c(3, 7, 12, 16))
res <- pileup(fl, scanBamParam=sbp, pileupParam=p_param)
xt <- xtabs(count ~ nucleotide + cycle_bin, res)
print(xt)
t(t(xt) / colSums(xt))

\donttest{
## single-width bins: count each position within a read separately.
p_param <- PileupParam(cycle_bins=seq(1,20))
res <- pileup(fl, scanBamParam=sbp, pileupParam=p_param)
xt <- xtabs(count ~ nucleotide + cycle_bin, res)
print(xt[ , c(1:3, 18:19)]) ## fit on one screen
print(t(t(xt) / colSums(xt))[ , c(1:3, 18:19)])
}
}