\name{FaFile}
\Rdversion{1.1}
\docType{class}
\alias{FaFile-class}
\alias{FaFileList-class}
% con/destructors
\alias{FaFile}
\alias{FaFileList}
\alias{open.FaFile}
\alias{close.FaFile }
% accessors
\alias{isOpen,FaFile-method}
% methods
\alias{indexFa,FaFile-method}
\alias{scanFaIndex,FaFile-method}
\alias{scanFaIndex,FaFileList-method}
\alias{countFa,FaFile-method}
\alias{scanFa,FaFile,GRanges-method}
\alias{scanFa,FaFile,RangesList-method}
\alias{scanFa,FaFile,RangedData-method}
\alias{scanFa,FaFile,missing-method}
\alias{getSeq,FaFile-method}
\alias{getSeq,FaFileList-method}

\title{Manipulate indexed fasta files.}

\description{

  Use \code{FaFile()} to create a reference to an indexed fasta
  file. The reference remains open across calls to methods, avoiding
  costly index re-loading.

  \code{FaFileList()} provides a convenient way of managing a list of
  \code{FaFile} instances.

}

\usage{

## Constructors

FaFile(file, ...)
FaFileList(...)

## Opening / closing

\S3method{open}{FaFile}(con, ...)
\S3method{close}{FaFile}(con, ...)

## accessors; also path(), index()

\S4method{isOpen}{FaFile}(con, rw="")

## actions

\S4method{indexFa}{FaFile}(file, ...)

\S4method{scanFaIndex}{FaFile}(file, ...)
\S4method{scanFaIndex}{FaFileList}(file, ..., as=c("GRangesList", "GRanges"))

\S4method{countFa}{FaFile}(file, ...)

\S4method{scanFa}{FaFile,GRanges}(file, param, ...)
\S4method{scanFa}{FaFile,RangesList}(file, param, ...)
\S4method{scanFa}{FaFile,RangedData}(file, param, ...)
\S4method{scanFa}{FaFile,missing}(file, param, ...)

\S4method{getSeq}{FaFile}(x, param, ...)
\S4method{getSeq}{FaFileList}(x, param, ...)
}

\arguments{

  \item{con, x}{An instance of \code{FaFile} or (for \code{getSeq})
    \code{FaFileList}.}

  \item{file}{A character(1) vector of the fasta file path (for
    \code{FaFile}), or an instance of class \code{FaFile} or
    \code{FaFileList} (for \code{scanFaIndex}, \code{getSeq}).}

  \item{param}{An optional \code{\linkS4class{GRanges}},
    \code{\linkS4class{RangesList}}, or \code{\linkS4class{RangedData}}
    instance to select reads (and sub-sequences) for input. See Methods,
    below.}

  \item{...}{Additional arguments. For \code{FaFileList}, this can
    either be a single character vector of paths to BAM files, or
    several instances of \code{FaFile} objects.}

  \item{rw}{Mode of file; ignored.}

  \item{as}{character(1) specifying the return type, selected from
  specified options. When \code{GRangesList}, index information from
  each file is returned as an element of the list. When
  \code{GRangesList}, index information is collapsed across files into
  the unique index elements.}

}

\section{Objects from the Class}{

  Objects are created by calls of the form \code{FaFile()}.

}

\section{Fields}{

  The \code{FaFile} class inherits fields from the
  \code{\linkS4class{RsamtoolsFile}} class.

}

\section{Functions and methods}{

  \code{FaFileList} inherits methods from
  \code{\link{RsamtoolsFileList}} and \code{\link{SimpleList}}.

  Opening / closing:
  \describe{

    \item{open.FaFile}{Opens the (local or remote) \code{path} and
      \code{index} files.  Returns a \code{FaFile} instance.}

    \item{close.FaFile}{Closes the \code{FaFile} \code{con}; returning
      (invisibly) the updated \code{FaFile}. The instance may be
      re-opened with \code{open.FaFile}.}

  }

  Accessors: 
  \describe{

    \item{path}{Returns a character(1) vector of the fasta path name.}

    \item{index}{Returns a character(1) vector of fasta index name
    (minus the '.fai' extension).}

  }

  Methods:
  \describe{

    \item{indexFa}{Visit the path in \code{path(file)} and create an
      index file (with the extension \sQuote{.fai}).}

    \item{scanFaIndex}{Read the sequence names and and widths of
      recorded in an indexed fasta file, returning the information as a
      \code{\linkS4class{GRanges}} object.}

    \item{countFa}{Return the number of records in the fasta file.}

    \item{scanFa}{Return the sequences indicated by \code{param} as a
      \code{\linkS4class{DNAStringSet}} instance. \code{seqnames(param)}
      selects the sequences to return; \code{start(param)} and
      \code{end{param}} define the (1-based) region of the sequence to
      return. Values of \code{end(param)} greater than the width of the
      sequence are set to the width of the sequence. When \code{param}
      is missing, all records are selected. When \code{length(param)==0}
      no records are selected.}

    \item{getSeq}{Returns the sequences indicated by \code{param} from
      the indexed fasta file(s) of \code{file}.

      For the \code{FaFile} method, the return type is a
      \code{DNAStringSet}. The \code{getSeq,FaFile} and
      \code{scanFa,FaFile,GRanges} methods differ in that \code{getSeq}
      will reverse complement sequences selected from the minus strand.

      For the \code{FaFileList} method, the \code{param} argument must
      be a \code{GRangesList} of the same length as \code{file},
      creating a one-to-one mapping between the ith element of
      \code{file} and the ith element of \code{param}; the return type
      is a \code{SimpleList} of \code{DNAStringSet} instances, with
      elements of the list in the same order as the input elements.}

    \item{show}{Compactly display the object.}

  }

}
\author{Martin Morgan}

\examples{

fl <- system.file("extdata", "ce2dict1.fa", package="Rsamtools")
fa <- open(FaFile(fl))                   # open
countFa(fa)
(idx <- scanFaIndex(fa))
(dna <- scanFa(fa, param=idx[1:2]))
ranges(idx) <- narrow(ranges(idx), -10)  # last 10 nucleotides
(dna <- scanFa(fa, param=idx[1:2]))

}

\keyword{classes}
