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
    ignore_query_Ns=TRUE, include_deletions=TRUE,
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

  The extended \code{CIGAR} format is partially supported by
  \code{pileup}, with a focus of supporting the range of deletion-like
  operations. The extended CIGAR allows a number of operations
  conceptually similar to a \sQuote{deletion} with respect to the
  reference genome, but offer more specialized meanings than a simple
  deletion. CIGAR \sQuote{N} operations (not to be confused with
  \sQuote{N} used for non-determinate bases in the \code{SEQ} field)
  indicate a large skipped region, \sQuote{S} a soft clip, and
  \sQuote{H} a hard clip. \sQuote{N}, \sQuote{S}, and \sQuote{H} CIGAR
  operations are never included in pileup counts. Only true deletions
  (\sQuote{D} in the CIGAR) can be included by setting
  \code{include_deletions} to \code{TRUE}.

  Insertion operations with respect to the reference genome (and
  therefore \sQuote{P} CIGAR operations) are ignored for \code{pileup}.

  \code{cycle_bins} allows users to differentiate pileup counts based on
  arbitrary regions within a read. \code{pileup} relies on
  \code{\link{cut}} to derive bins, but limits input to non-negative
  numeric values (coerced to \code{integer}s), including Inf. If
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

}

\value{
  
  \code{pileup} returns a \code{data.frame} with columns
  \code{seqnames}, \code{pos}, and \code{count} minimally, and optional
  \code{strand}, \code{nulceotide}, and \code{cycle_bin} columns as
  dictated by arguments to \code{PileupParam}.

  If a \code{which} argument is specified for the scanBamParam, a
  \code{which_label} column (in the form \sQuote{rname:start-end}) is
  automatically included. The \code{which_label} column is used to
  maintain grouping of results, such that two queries of the same
  genomic region can be differentiated.

  Order of rows in \code{data.frame} is first by order of
  \code{seqnames} column based on the BAM index file, then
  non-decreasing order on columns \code{pos}, then \code{nucleotide},
  then \code{strand}.

  \code{PileupParam} returns an instance of PileupParam class.

}

\seealso{
  \itemize{
    \item{\link{Rsamtools}}
    \item{\link{ScanBamParam}}
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

## pileup for all of chr14 with default Pileup parameters
library("RNAseqData.HNRNPC.bam.chr14")

fl <- RNAseqData.HNRNPC.bam.chr14_BAMFILES[1]
sbp <- ScanBamParam(which=GRanges("chr14", IRanges(1, 107349540)))
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

\dontrun{
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

\dontrun{
## single-width bins: count each position within a read separately.
p_param <- PileupParam(cycle_bins=seq(1,20))
res <- pileup(fl, scanBamParam=sbp, pileupParam=p_param)
xt <- xtabs(count ~ nucleotide + cycle_bin, res)
head(print(xt))
head(t(t(xt) / colSums(xt)))
}
}