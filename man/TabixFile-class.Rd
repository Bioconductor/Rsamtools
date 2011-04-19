\name{TabixFile}
\Rdversion{1.1}
\docType{class}
\alias{TabixFile-class}
% con/destructors
\alias{TabixFile}
\alias{open.TabixFile}
\alias{close.TabixFile }
% accessors
\alias{isOpen,TabixFile-method}
% methods
\alias{seqnamesTabix,TabixFile-method}
\alias{scanTabix,TabixFile,RangesList-method}
\alias{scanTabix,TabixFile,RangedData-method}
\alias{scanTabix,TabixFile,GRanges-method}
\alias{yieldTabix,TabixFile-method}

\title{Manipulate tabix indexed tab-delimited files.}

\description{

  Use \code{TabixFile()} to create a reference to a Tabix file (and its
  index). Once opened, the reference remains open across calls to
  methods, avoiding costly index re-loading.

}

\usage{

## Opening / closing

TabixFile(file, index = paste(file, "tbi", sep="."), ...)
\S3method{open}{TabixFile}(con, ...)
\S3method{close}{TabixFile}(con, ...)

## accessors; also path(), index()

\S4method{isOpen}{TabixFile}(con, rw="")

## actions

\S4method{seqnamesTabix}{TabixFile}(file, ...)
\S4method{scanTabix}{TabixFile,RangesList}(file, ..., param)
\S4method{scanTabix}{TabixFile,RangedData}(file, ..., param)
\S4method{scanTabix}{TabixFile,GRanges}(file, ..., param)

\S4method{yieldTabix}{TabixFile}(file, ..., yieldSize=1000000L)

}

\arguments{

  \item{con}{An instance of \code{TabixFile}.}

  \item{file}{For TabixFile(), A character(1) vector to the tabix file
    path; can be remote (http://, ftp://). For others, a
    \code{TabixFile} instance.}

  \item{index}{A character(1) vector of the tabix file index.}

  \item{param}{An instance of GRanges, IRangedData, or RangesList, used
    to select which records to scan.}

  \item{yieldSize}{integer(1) indicating the maximum number of records
    to retreieve.}

  \item{...}{Additional arguments, currently unused.}

  \item{rw}{character() indicating mode of file; not used for \code{TabixFile}.}

}

\section{Objects from the Class}{

  Objects are created by calls of the form \code{TabixFile()}.

}

\section{Fields}{

  The \code{TabixFile} class inherits fields from the
  \code{\linkS4class{RsamtoolsFile}} class.

}

\section{Functions and methods}{

  Opening / closing:
  \describe{

    \item{open.TabixFile}{Opens the (local or remote) \code{path} and
      \code{index}. Returns a \code{TabixFile} instance.}

    \item{close.TabixFile}{Closes the \code{TabixFile} \code{con}; returning
	  (invisibly) the updated \code{TabixFile}. The instance may be
	  re-opened with \code{open.TabixFile}.}

  }

  Accessors: 
  \describe{

    \item{path}{Returns a character(1) vector of the tabix path name.}

    \item{index}{Returns a character(1) vector of tabix index name.}

  }

  Methods:
  \describe{

    \item{seqnamesTabix}{Visit the path in \code{path(file)}, returning
      the sequence names present in the file.}

    \item{scanTabix}{Visit the path in \code{path(file)}, returning the
      result of \code{\link{scanTabix}} applied to the specified path.}

    \item{show}{Compactly display the object.}

  }

}
\author{Martin Morgan}

\examples{

fl <- system.file("extdata", "example.gtf.gz", package="Rsamtools")
tbx <- TabixFile(fl)

param <- GRanges(c("chr1", "chr2"), IRanges(c(1, 1), width=100000))
res <- scanTabix(tbx, param=param)
names(res)
res[["chr1:1-100000"]][1:2]

## parse 100 records at a time
tbx <- open(TabixFile(fl))
while(length(res <- yieldTabix(tbx, yieldSize=100L)))
   cat("records read:", length(res), "\n")
close(tbx)

}

\keyword{classes}
