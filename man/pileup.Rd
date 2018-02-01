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
\alias{left_bins}
\alias{query_bins}
\alias{cycle_bins}

% pileup
\alias{pileup}
\alias{pileup,character-method}
\alias{pileup,BamFile-method}

% utility
\alias{phred2ASCIIOffset}

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
PileupParam(max_depth=250, min_base_quality=13, min_mapq=0,
    min_nucleotide_depth=1, min_minor_allele_depth=0,
    distinguish_strands=TRUE, distinguish_nucleotides=TRUE,
    ignore_query_Ns=TRUE, include_deletions=TRUE, include_insertions=FALSE,
    left_bins=NULL, query_bins=NULL, cycle_bins=NULL)

phred2ASCIIOffset(phred=integer(),
    scheme= c("Illumina 1.8+", "Sanger", "Solexa", "Illumina 1.3+",
              "Illumina 1.5+"))
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
    each nucleotide in an alignment. Use \code{phred2ASCIIOffset} to
    help translate numeric or character values to these offsets.}

  \item{min_mapq}{integer(1); minimum \sQuote{MAPQ} value for an
    alignment to be included in pileup.}

  \item{min_nucleotide_depth}{integer(1); minimum count of each
    nucleotide (\emph{independent} of other nucleotides) at a given
    position required for said nucleotide to appear in the result.}

  \item{min_minor_allele_depth}{integer(1); minimum count of \emph{all}
    nucleotides other than the major allele at a given position required
    for a particular nucleotide to appear in the result.}

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

  \item{left_bins}{numeric; all same sign; unique positions within a
    read to delimit bins. Position within read is based on counting from
    the \emph{5' end regardless of strand}. Minimum of two values are
    required so at least one range can be formed. \code{NULL} (default)
    indicates no binning. Use negative values to count from the opposite
    end. Sorted order not required. Floating-point values are coerced to
    \code{integer}.

    If you want to delimit bins based on sequencing cycles to, e.g.,
    discard later cycles, \code{\link{query_bins}} probably gives the
    desired behavior.

  }

  \item{query_bins}{numeric; all same sign; unique positions within a
    read to delimit bins. Position within a read is based on counting
    from the \emph{5' end} when the read aligns to \emph{plus strand},
    counting from the \emph{3' end} when read aligns to minus
    strand. Minimum of two values are required so at least one range can
    be formed. \code{NULL} (default) indicates no binning. Use negative
    values to count from the opposite end. Sorted order not
    required. Floating-point values are coerced to \code{integer}.}

  \item{phred}{Either a numeric() (coerced to integer()) PHRED score
    (e.g., in the range (0, 41) for the \sQuote{Illumina 1.8+} scheme)
    or character() of printable ASCII characters. When \code{phred} is
    character(), it can be of length 1 with 1 or more characters, or of
    any length where all elements have exactly 1 character.}

  \item{scheme}{Encoding scheme, used to translate numeric() PHRED
    scores to their ASCII code. Ignored when \code{phred} is already
    character().}

  \item{cycle_bins}{DEPRECATED. See \code{\link{left_bins}} for
    identical behavior.}

}

\details{

  \emph{NB}: Use of \code{pileup} assumes familiarity with
  \code{\link{ScanBamParam}}, and use of \code{left_bins} and
  \code{query_bins} assumes familiarity with \code{\link{cut}}.

  \code{pileup} visits each position in the BAM file, subject to
  constraints implied by \code{which} and \code{flag} of
  \code{scanBamParam}. For a given position, all reads overlapping the
  position that are consistent with constraints dictated by \code{flag}
  of \code{scanBamParam} and \code{pileupParam} are used for
  counting. \code{pileup} also automatically excludes reads with flags
  that indicate any of:

  \itemize{
    \item{unmapped alignment (\code{\link{isUnmappedQuery}})}
    \item{secondary alignment (\code{\link{isSecondaryAlignment}})}
    \item{not passing quality controls (\code{\link{isNotPassingQualityControls}})}
    \item{PCR or optical duplicate (\code{\link{isDuplicate}})}
   }

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
      \sQuote{H} a hard clip. \sQuote{N}, \sQuote{S}, and \sQuote{H}
      CIGAR operations are never counted: only counts of true deletions
      (\sQuote{D} in the CIGAR) can be included by setting
      \code{include_deletions} to \code{TRUE}.

      Soft clip operations contribute to the relative position of
      nucleotides within a read, whereas hard clip and refskip
      operations do not. For example, consider a sequence in a bam file
      that starts at 11, with a CIGAR \code{2S1M} and \code{SEQ} field
      \code{ttA}. The cycle position of the \code{A} nucleotide will be
      \code{3}, but the reported position will be \code{11}. If using
      \code{left_bins} or \code{query_bins} it might make sense to first
      preprocess your files to remove soft clipping.

    }

    \item{Insertions and padding:

      \code{pileup} can include counts of insertion operations by
      setting \code{include_insertions} to \code{TRUE}. Details about
      insertions are effectively truncated such that each insertion is
      reduced to a single \sQuote{+} nucleotide. Information about
      e.g. the nucleotide code and base quality within the insertion is
      not included.

      Because \sQuote{+} is used as a nucleotide code to represent
      insertions in \code{pileup}, counts of the \sQuote{+} nucleotide
      participate in voting for \code{min_nucleotide_depth} and
      \code{min_minor_allele_depth} functionality.
      
      The position of an insertion is the position that precedes it on
      the reference sequence. \emph{Note:} this means if
      \code{include_insertions} is \code{TRUE} a single read will
      contribute two nucleotide codes (\code{+}, plus that of the
      non-insertion base) at a given position if the non-insertion base
      passes filter criteria.

      \sQuote{P} CIGAR (padding) operations never affect counts.

    }
  }

  The \dQuote{bin} arguments \code{query_bins} and \code{left_bins}
  allow users to differentiate pileup counts based on arbitrary
  non-overlapping regions within a read. \code{pileup} relies on
  \code{\link{cut}} to derive bins, but limits input to numeric values
  of the same sign (coerced to \code{integer}s), including
  \code{+/-Inf}. If a \dQuote{bin} argument is set \code{pileup}
  automatically excludes bases outside the implicit outer range. Here
  are some important points regarding bin arguments:

  \itemize{

    \item{\code{query_bins} vs. \code{left_bins}:}{

      BAM files store sequence data from the 5' end to the 3' end
      (regardless of to which strand the read aligns). That means for
      reads aligned to the minus strand, cycle position within a read is
      effectively reversed. Take care to use the appropriate bin
      argument for your use case.

      The most common use of a bin argument is to bin later cycles
      separately from earlier cycles; this is because accuracy typically
      degrades in later cycles. For this application, \code{query_bins}
      yields the correct behavior because bin positions are relative to
      cycle order (i.e., sensitive to strand).

      \code{left_bins} (in contrast with \code{query_bins}) determines
      bin positions from the 5' end, regardless of strand.

    }

    \item{Positive or negative bin values can be used to delmit bins
      based on the effective \dQuote{start} or \dQuote{end} of a
      read. For \code{left_bin} the effective start is always the 5' end
      (left-to-right as appear in the BAM file).

      For \code{query_bin} the effective start is the first cycle of the
      read as it came off the sequencer; that means the 5' end for reads
      aligned to the plus strand and 3' for reads aligned to the minus
      strand.

      \itemize{

        \item{\emph{From effective start of reads}: use positive values,
          \code{0}, and \code{(+)Inf}. Because \code{\link{cut}}
          produces ranges in the form (first,last], \sQuote{0} should be
          used to create a bin that includes the first position. To
          account for variable-length reads in BAM files, use
          \code{(+)Inf} as the upper bound on a bin that extends to the
          end of all reads.}

        \item{\emph{From effective end of reads}: use negative values
          and \code{-Inf}. \code{-1} denotes the last position in a
          read. Because \code{\link{cut}} produces ranges in the form
          (first,last], specify the lower bound of a bin by using one
          less than the desired value, e.g., a bin that captures only
          the second nucleotide from the last position:
          \code{query_bins=c(-3, -2)}. To account for variable-length
          reads in BAM files, use \code{-Inf} as the lower bound on a
          bin that extends to the beginning of all reads.}
      }
    }
  } %% end subsection bin argument gotchas

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
  \code{strand}, \code{nulceotide}, and \code{left_bin} /
  \code{query_bin} columns appear as dictated by arguments to
  \code{PileupParam}.

  Column details:

  \itemize{
    
    \item{\code{seqnames}: factor. SAM \sQuote{RNAME} field of
    alignment.}

    \item{\code{pos}: integer(1). Genomic position of base. Derived by
    offset from SAM \sQuote{POS} field of alignment.}

    \item{\code{strand}: factor. \sQuote{strand} to which read aligns.}

    \item{\code{nucleotide}: factor. \sQuote{nucleotide} of base,
      extracted from SAM \sQuote{SEQ} field of alignment.}

    \item{\code{left_bin} / \code{query_bin}: factor. Bin in which base
      appears.}

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
  then \code{strand}, then \code{left_bin} / \code{query_bin}.

  \code{PileupParam} returns an instance of PileupParam class.

  \code{phred2ASCIIOffset} returns a named integer vector of the same
  length or number of characters in \code{phred}. The names are the
  ASCII encoding, and the values are the offsets to be used with
  \code{min_base_quality}.

}

\seealso{
  \itemize{
    \item{\link{Rsamtools}}
    \item{\linkS4class{BamFile}}
    \item{\link{ScanBamParam}}
    \item{\link{scanBam}}
    \item{\link{cut}}
  }

  For the relationship between PHRED scores and ASCII encoding, see
  \url{https://en.wikipedia.org/wiki/FASTQ_format#Encoding}.
}

\author{Nathaniel Hayden \url{nhayden@fredhutch.org}}

%% This information should go somewhere--essential for an outside
%% developer to help build context quickly

%% \excruciatingdetails{

%%   There are really three stages of filtering; (1) filtering based on
%%   \code{scanBamParam} arguments (via _do_scan_bam, (2) filtering
%%   based on arguments of \code{pileupParam} that apply on a per-read
%%   and per-alignment basis and lastly, (3) a final filtering stage
%%   applied on the intermediate count of a position after all
%%   overlapping alignments have been
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
## Minimum base quality to be tallied
p_param <- PileupParam(min_base_quality = 10L)

\donttest{
## no 'which' argument to ScanBamParam: process entire file at once
res <- pileup(fl, pileupParam = p_param)
head(res)
table(res$strand, res$nucleotide)
}

\donttest{
## no 'which' argument to ScanBamParam with 'yieldSize' set on BamFile
bf <- open(BamFile(fl, yieldSize=5e4))
repeat {
    res <- pileup(bf, pileupParam = p_param)
    message(nrow(res), " rows in result data.frame")
    if(nrow(res) == 0L)
        break
}
close(bf)
}

## pileup for the first half of chr14 with default Pileup parameters
## 'which' argument: process only specified genomic range(s)
sbp <- ScanBamParam(which=GRanges("chr14", IRanges(1, 53674770)))
res <- pileup(fl, scanBamParam=sbp, pileupParam = p_param)
head(res)
table(res$strand, res$nucleotide)

## specify pileup parameters: include ambiguious nucleotides
## (the 'N' level in the nucleotide column of the data.frame)
p_param <- PileupParam(ignore_query_Ns=FALSE, min_base_quality = 10L)
res <- pileup(fl, scanBamParam=sbp, pileupParam=p_param)
head(res)
table(res$strand, res$nucleotide)

## Don't distinguish strand, require a minimum frequency of 7 for a
## nucleotide at a genomic position to be included in results.

p_param <- PileupParam(distinguish_strands=TRUE,
                       min_base_quality = 10L,
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

## ASCII offsets to help specify min_base_quality, e.g., quality of at
## least 10 on the Illumina 1.3+ scale
phred2ASCIIOffset(10, "Illumina 1.3+")

## Well-supported polymorphic positions (minor allele present in at
## least 5 alignments) with high map quality
p_param <- PileupParam(min_minor_allele_depth=5,
                       min_mapq=40,
                       min_base_quality = 10,
                       distinguish_strand=FALSE)
res <- pileup(fl, scanBamParam=sbp, pileupParam=p_param)
dim(res) ## reduced to our biologically interesting positions
head(xtabs(count ~ pos + nucleotide, res))

## query_bins

\donttest{

## basic use of positive bins: single pivot; count bases that appear in
## first 10 cycles of a read separately from the rest
p_param <- PileupParam(query_bins=c(0, 10, Inf), min_base_quality = 10)
res <- pileup(fl, scanBamParam=sbp, pileupParam=p_param)

## basic use of positive bins: simple range; only include bases in
## cycle positions 6-10 within read
p_param <- PileupParam(query_bins=c(5, 10), min_base_quality = 10)
res <- pileup(fl, scanBamParam=sbp, pileupParam=p_param)

## basic use of negative bins: single pivot; count bases that appear in
## last 3 cycle positions of a read separately from the rest.
p_param <- PileupParam(query_bins=c(-Inf, -4, -1), min_base_quality = 10)
res <- pileup(fl, scanBamParam=sbp, pileupParam=p_param)

## basic use of negative bins: drop nucleotides in last two cycle
## positions of each read
p_param <- PileupParam(query_bins=c(-Inf, -3), min_base_quality = 10)
res <- pileup(fl, scanBamParam=sbp, pileupParam=p_param)
}

## typical use: beginning, middle, and end segments; because of the
## nature of sequencing technology, it is common for bases in the
## beginning and end segments of each read to be less reliable. pileup
## makes it easy to count each segment separately.

## Assume determined ahead of time that for the data 1-7 makes sense for
## beginning, 8-12 middle, >=13 end (actual cycle positions should be
## tailored to data in actual BAM files).
p_param <- PileupParam(query_bins=c(0, 7, 12, Inf), ## note non-symmetric
                       min_base_quality = 10)
res <- pileup(fl, scanBamParam=sbp, pileupParam=p_param)
xt <- xtabs(count ~ nucleotide + query_bin, res)
print(xt)
t(t(xt) / colSums(xt)) ## cheap normalization for illustrative purposes

## 'implicit outer range': in contrast to c(0, 7, 12, Inf), suppose we
##  still want to have beginning, middle, and end segements, but know
##  that the first three cycles and any bases beyond the 16th cycle are
##  irrelevant. Hence, the implicit outer range is (3,16]; all bases
##  outside of that are dropped.
p_param <- PileupParam(query_bins=c(3, 7, 12, 16), min_base_quality = 10)
res <- pileup(fl, scanBamParam=sbp, pileupParam=p_param)
xt <- xtabs(count ~ nucleotide + query_bin, res)
print(xt)
t(t(xt) / colSums(xt))

\donttest{
## single-width bins: count each cycle within a read separately.
p_param <- PileupParam(query_bins=seq(1,20), min_base_quality = 10)
res <- pileup(fl, scanBamParam=sbp, pileupParam=p_param)
xt <- xtabs(count ~ nucleotide + query_bin, res)
print(xt[ , c(1:3, 18:19)]) ## fit on one screen
print(t(t(xt) / colSums(xt))[ , c(1:3, 18:19)])
}
}