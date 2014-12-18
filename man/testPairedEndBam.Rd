\name{testPairedEndBam}
\alias{testPairedEndBam}
\alias{testPairedEndBam,character-method}
\alias{testPairedEndBam,BamFile-method}

\title{
    Quickly test if a BAM file has paired end reads
}

\description{

    Iterate through a BAM file until a paired-end read is encountered or
    the end of file is reached; report the occurrence of paired-end
    reads to the user.

}

\usage{
testPairedEndBam(file, index=file, ...)
}

\arguments{

    \item{file}{character(1) BAM file name, or a \code{\link{BamFile}}
      instance. Open \code{BamFile}s are closed; their yield size is
      respected when iterating through the file.}
        
    \item{index}{(optional) character(1) name of the index file of the
        'BAM' file being processed; this is given \emph{without} the
        '.bai' extension.}
    
    \item{\dots}{Additional arguments, currently unused.}

}

\value{

    A logical vector of length 1 containing TRUE is returned if BAM file 
    contained paired end reads, FALSE otherwise.  

}

\author{

  Martin Morgan \url{mailto:mtmorgan@fhcrc.org}, Sonali Arora
  \url{mailto:sarora@fhcrc.org}

}

\examples{
fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
testPairedEndBam(fl)
}

