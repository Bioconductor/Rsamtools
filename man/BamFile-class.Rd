\name{BamFile}
\Rdversion{1.1}
\docType{class}
\alias{BamFile-class}
\alias{BamFileList-class}
% con/destructors
\alias{BamFile}
\alias{BamFileList}
\alias{open.BamFile}
\alias{close.BamFile}
% accessors
\alias{isOpen,BamFile-method}
\alias{isIncomplete,BamFile-method}
% methods
\alias{scanBamHeader,BamFile-method}
\alias{seqinfo,BamFile-method}
\alias{seqinfo,BamFileList-method}
\alias{obeyQname}
\alias{obeyQname<-}
\alias{obeyQname,BamFile-method}
\alias{obeyQname<-,BamFile-method}
\alias{obeyQname,BamFileList-method}
\alias{obeyQname<-,BamFileList-method}
\alias{asMates}
\alias{asMates<-}
\alias{asMates,BamFile-method}
\alias{asMates<-,BamFile-method}
\alias{asMates,BamFileList-method}
\alias{asMates<-,BamFileList-method}
\alias{qnamePrefixEnd}
\alias{qnamePrefixEnd<-}
\alias{qnamePrefixEnd,BamFile-method}
\alias{qnamePrefixEnd<-,BamFile-method}
\alias{qnamePrefixEnd,BamFileList-method}
\alias{qnamePrefixEnd<-,BamFileList-method}
\alias{qnameSuffixStart}
\alias{qnameSuffixStart<-}
\alias{qnameSuffixStart,BamFile-method}
\alias{qnameSuffixStart<-,BamFile-method}
\alias{qnameSuffixStart,BamFileList-method}
\alias{qnameSuffixStart<-,BamFileList-method}
\alias{scanBam,BamFile-method}
\alias{idxstatsBam,BamFile-method}
\alias{countBam,BamFile-method}
\alias{countBam,BamFileList-method}
\alias{filterBam,BamFile-method}
\alias{indexBam,BamFile-method}
\alias{sortBam,BamFile-method}

\alias{mergeBam,BamFileList-method}
% counting
\alias{quickBamFlagSummary,BamFile-method}

\alias{show,BamFile-method}
\alias{show,BamFileList-method}

\title{Maintain and use BAM files}

\description{

  Use \code{BamFile()} to create a reference to a BAM file (and
  optionally its index). The reference remains open across calls to
  methods, avoiding costly index re-loading.

  \code{BamFileList()} provides a convenient way of managing a list of
  \code{BamFile} instances.

}

\usage{

## Constructors

BamFile(file, index=file, ..., yieldSize=NA_integer_, obeyQname=FALSE,
        asMates=FALSE, qnamePrefixEnd=NA, qnameSuffixStart=NA)
BamFileList(..., yieldSize=NA_integer_, obeyQname=FALSE, asMates=FALSE,
            qnamePrefixEnd=NA, qnameSuffixStart=NA)

## Opening / closing

\S3method{open}{BamFile}(con, ...)
\S3method{close}{BamFile}(con, ...)

## accessors; also path(), index(), yieldSize()

\S4method{isOpen}{BamFile}(con, rw="")
\S4method{isIncomplete}{BamFile}(con)
\S4method{obeyQname}{BamFile}(object, ...)
obeyQname(object, ...) <- value
\S4method{asMates}{BamFile}(object, ...)
asMates(object, ...) <- value
\S4method{qnamePrefixEnd}{BamFile}(object, ...)
qnamePrefixEnd(object, ...) <- value
\S4method{qnameSuffixStart}{BamFile}(object, ...)
qnameSuffixStart(object, ...) <- value

## actions

\S4method{scanBamHeader}{BamFile}(files, ..., what=c("targets", "text"))
\S4method{seqinfo}{BamFile}(x)
\S4method{seqinfo}{BamFileList}(x)
\S4method{filterBam}{BamFile}(file, destination, index=file, ...,
    filter=FilterRules(), indexDestination=TRUE,
    param=ScanBamParam(what=scanBamWhat()))
\S4method{indexBam}{BamFile}(files, ...)
\S4method{sortBam}{BamFile}(file, destination, ..., byQname=FALSE, maxMemory=512, byTag=NULL, nThreads=1L)
\S4method{mergeBam}{BamFileList}(files, destination, ...)

## reading

\S4method{scanBam}{BamFile}(file, index=file, ..., param=ScanBamParam(what=scanBamWhat()))

## counting

\S4method{idxstatsBam}{BamFile}(file, index=file, ...)
\S4method{countBam}{BamFile}(file, index=file, ..., param=ScanBamParam())
\S4method{countBam}{BamFileList}(file, index=file, ..., param=ScanBamParam())
\S4method{quickBamFlagSummary}{BamFile}(file, ..., param=ScanBamParam(), main.groups.only=FALSE)
}

\arguments{

    \item{...}{Additional arguments.

      For \code{BamFileList}, this can either be a single character vector
      of paths to BAM files, or several instances of \code{BamFile}
      objects. When a character vector of paths, a second named argument
      \sQuote{index} can be a \code{character()} vector of length equal to
      the first argument specifying the paths to the index files, or
      \code{character()} to indicate that no index file is available. See
      \code{\link{BamFile}}.
    }

    \item{con}{An instance of \code{BamFile}.}

    \item{x, object, file, files}{A character vector of BAM file paths
      (for \code{BamFile}) or a \code{BamFile} instance (for other
      methods).}

    \item{index}{character(1); the BAM index file path (for
      \code{BamFile}); ignored for all other methods on this page.}

    \item{yieldSize}{Number of records to yield each time the file
      is read from with \code{scanBam}. See \sQuote{Fields}
      section for details.}

    \item{asMates}{Logical indicating if records should be paired
      as mates. See \sQuote{Fields} section for details.}

    \item{qnamePrefixEnd}{Single character (or NA) marking the 
      end of the qname prefix. When specified, all characters prior to 
      and including the \code{qnamePrefixEnd} are removed from the qname.
      If the prefix is not found in the qname the qname is not trimmed.
      Currently only implemented for mate-pairing (i.e., when
      \code{asMates=TRUE} in a BamFile.}

    \item{qnameSuffixStart}{Single character (or NA) marking the 
      start of the qname suffix. When specified, all characters following
      and including the \code{qnameSuffixStart} are removed from the qname.
      If the suffix is not found in the qname the qname is not trimmmed.
      Currently only implemented for mate-pairing (i.e., when
      \code{asMates=TRUE} in a BamFile.}

    \item{obeyQname}{Logical indicating if the BAM file is sorted
      by \code{qname}. In Bioconductor > 2.12 paired-end files do
      not need to be sorted by \code{qname}. Instead use
      \code{asMates=TRUE} for reading paired-end data. See 
      \sQuote{Fields} section for details.}

    \item{value}{Logical value for setting \code{asMates} and
      \code{obeyQname} in a BamFile instance.}

    \item{what}{For \code{scanBamHeader}, a character vector specifying
      that either or both of \code{c("targets", "text")} are to be
      extracted from the header; see \code{\link{scanBam}} for
      additional detail.}

    \item{filter}{A \code{\link{FilterRules}} instance. Functions in the
      \code{FilterRules} instance should expect a single \code{DataFrame}
      argument representing all information specified by
      \code{param}. Each function must return a \code{logical} vector,
      usually of length equal to the number of rows of the
      \code{DataFrame}. Return values are used to include (when
      \code{TRUE}) corresponding records in the filtered BAM file.}

    \item{destination}{character(1) file path to write filtered reads to.}

    \item{indexDestination}{logical(1) indicating whether the destination
      file should also be indexed.}

    \item{byQname, maxMemory, byTag, nThreads}{See \code{\link{sortBam}}.}

    \item{param}{An optional \code{\linkS4class{ScanBamParam}} instance to
       further influence scanning, counting, or filtering.}

    \item{rw}{Mode of file; ignored.}

    \item{main.groups.only}{See \code{\link{quickBamFlagSummary}}.}
}

\section{Objects from the Class}{

  Objects are created by calls of the form \code{BamFile()}.

}

\section{Fields}{

  The \code{BamFile} class inherits fields from the
  \code{\linkS4class{RsamtoolsFile}} class and has fields:

  \describe{

    \item{yieldSize: }{Number of records to yield each time the file is
      read from using \code{scanBam} or, when \code{length(bamWhich())
      != 0}, a threshold which yields records in complete ranges whose
      sum first exceeds \code{yieldSize}. Setting \code{yieldSize} on a
      \code{BamFileList} does not alter existing yield sizes set on the
      individual \code{BamFile} instances.}

    \item{asMates: }{A logical indicating if the records should be
      returned as mated pairs. When \code{TRUE} \code{scanBam} attempts
      to mate (pair) the records and returns two additional fields
      \code{groupid} and \code{mate_status}. \code{groupid} is an integer
      vector of unique group ids; \code{mate_status} is a factor with level
      \code{mated} for records successfully paired by the algorithm,
      \code{ambiguous} for records that are possibly mates but cannot be
      assigned unambiguously, or \code{unmated} for reads that did not
      have valid mates.

      Mate criteria:
      \itemize{
        \item Bit 0x40 and 0x80: Segments are a pair of first/last OR
              neither segment is marked first/last
        \item Bit 0x100: Both segments are secondary OR both not secondary
        \item Bit 0x10 and 0x20: Segments are on opposite strands
        \item mpos match:
              segment1 mpos matches segment2 pos AND
              segment2 mpos matches segment1 pos
        \item tid match
      }

      Flags, tags and ranges may be specified in the \code{ScanBamParam}
      for fine tuning of results.}

      \item{obeyQname: }{A logical(0) indicating if the file was sorted by 
        qname. In Bioconductor > 2.12 paired-end files do not need to be 
        sorted by \code{qname}. Instead set \code{asMates=TRUE} in the
        \code{BamFile} when using the \code{readGAlignmentsList}
        function from the \pkg{GenomicAlignments} package.
      }
  }
}

\section{Functions and methods}{

  \code{BamFileList} inherits additional methods from
  \code{\link{RsamtoolsFileList}} and \code{\link{SimpleList}}.

  Opening / closing:
  \describe{

    \item{open.BamFile}{Opens the (local or remote) \code{path} and
      \code{index} (if \code{bamIndex} is not \code{character(0)}),
      files.  Returns a \code{BamFile} instance.}

    \item{close.BamFile}{Closes the \code{BamFile} \code{con}; returning
      (invisibly) the updated \code{BamFile}. The instance may be
      re-opened with \code{open.BamFile}.}

    \item{isOpen}{Tests whether the \code{BamFile} \code{con} has been
      opened for reading.}

    \item{isIncomplete}{Tests whether the \code{BamFile} \code{con} is
      niether closed nor at the end of the file.}

  }

  Accessors: 
  \describe{

    \item{path}{Returns a character(1) vector of BAM path names.}

    \item{index}{Returns a character(0) or character(1) vector of BAM
      index path names.}

    \item{yieldSize, yieldSize<-}{Return or set an integer(1) vector
      indicating yield size.}

    \item{obeyQname, obeyQname<-}{Return or set a logical(0)
      indicating if the file was sorted by qname.}

    \item{asMates, asMates<-}{Return or set a logical(0)
      indicating if the records should be returned as mated pairs.}
  }

  Methods:
  \describe{

    \item{scanBamHeader}{Visit the path in \code{path(file)}, returning
      the information contained in the file header; see
      \code{\link{scanBamHeader}}.}

    \item{seqinfo, seqnames, seqlength}{Visit the path in
      \code{path(file)}, returning a \code{\linkS4class{Seqinfo}},
      character, or named integer vector containing information on the
      anmes and / or lengths of each sequence. Seqnames are ordered
      as they appear in the file.}

    \item{scanBam}{Visit the path in \code{path(file)}, returning the
      result of \code{\link{scanBam}} applied to the specified path.}

    \item{countBam}{Visit the path(s) in \code{path(file)}, returning
      the result of \code{\link{countBam}} applied to the specified
      path.}

    \item{idxstatsBam}{Visit the index in \code{index(file)}, quickly
      returning a \code{data.frame} with columns \code{seqnames},
      \code{seqlength}, \code{mapped} (number of mapped reads on
      \code{seqnames}) and \code{unmapped} (number of unmapped reads).}

    \item{filterBam}{Visit the path in \code{path(file)}, returning the
      result of \code{\link{filterBam}} applied to the specified
      path. A single file can be filtered to one or several
      destinations, as described in \code{\link{filterBam}}.}

    \item{indexBam}{Visit the path in \code{path(file)}, returning
      the result of \code{\link{indexBam}} applied to the specified
      path.}

    \item{sortBam}{Visit the path in \code{path(file)}, returning the
      result of \code{\link{sortBam}} applied to the specified path.}

    \item{mergeBam}{Merge several BAM files into a single BAM file. See
      \code{\link{mergeBam}} for details; additional arguments supported
      by \code{mergeBam,character-method} are also available for
      \code{BamFileList}.}

    \item{show}{Compactly display the object.}

  }

}

\author{Martin Morgan and Marc Carlson}

\seealso{
  \itemize{
    \item The \code{\link[GenomicAlignments]{readGAlignments}},
          \code{\link[GenomicAlignments]{readGAlignmentPairs}},
          and \code{\link[GenomicAlignments]{readGAlignmentsList}}
          functions defined in the \pkg{GenomicAlignments} package.

    \item \code{\link[GenomicAlignments]{summarizeOverlaps}} and
          \link[GenomicAlignments]{findSpliceOverlaps-methods} in the
          \pkg{GenomicAlignments} package for methods that work on a
          \link{BamFile} and \link{BamFileList} objects.
  }
}

\examples{

##
## BamFile options.
##

fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
bf <- BamFile(fl)
bf

## When 'asMates=TRUE' scanBam() reads the data in as
## pairs. See 'asMates' above for details of the pairing
## algorithm.
asMates(bf) <- TRUE

## When 'yieldSize' is set, scanBam() will iterate
## through the file in chunks.
yieldSize(bf) <- 500 

## Some applications append a filename (e.g., NCBI Sequence Read 
## Archive (SRA) toolkit) or allele identifier to the sequence qname.
## This may result in a unique qname for each record which presents a
## problem when mating paired-end reads (identical qnames is one
## criteria for paired-end mating). 'qnamePrefixEnd' and 
## 'qnameSuffixStart' can be used to trim an unwanted prefix or suffix.
qnamePrefixEnd(bf) <- "/"
qnameSuffixStart(bf) <- "." 

##
## Reading Bam files.
##

fl <- system.file("extdata", "ex1.bam", package="Rsamtools",
                  mustWork=TRUE)
(bf <- BamFile(fl))
head(seqlengths(bf))                    # sequences and lengths in BAM file

if (require(RNAseqData.HNRNPC.bam.chr14)) {
    bfl <- BamFileList(RNAseqData.HNRNPC.bam.chr14_BAMFILES)
    bfl
    bfl[1:2]                            # subset
    bfl[[1]]                            # select first element -- BamFile
    ## merged across BAM files
    seqinfo(bfl)
    head(seqlengths(bfl))
}


length(scanBam(fl)[[1]][[1]])  # all records

bf <- open(BamFile(fl))        # implicit index
bf
identical(scanBam(bf), scanBam(fl))
close(bf)

## Use 'yieldSize' to iterate through a file in chunks.
bf <- open(BamFile(fl, yieldSize=1000)) 
while (nrec <- length(scanBam(bf)[[1]][[1]]))
    cat("records:", nrec, "\n")
close(bf)

## Repeatedly visit multiple ranges in the BamFile. 
rng <- GRanges(c("seq1", "seq2"), IRanges(1, c(1575, 1584)))
bf <- open(BamFile(fl))
sapply(seq_len(length(rng)), function(i, bamFile, rng) {
    param <- ScanBamParam(which=rng[i], what="seq")
    bam <- scanBam(bamFile, param=param)[[1]]
    alphabetFrequency(bam[["seq"]], baseOnly=TRUE, collapse=TRUE)
}, bf, rng)
close(bf)

}

\keyword{classes}
