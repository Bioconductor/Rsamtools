\name{testPairedEndBam}
\alias{testPairedEndBam}
\alias{testPairedEndBam,character-method}
\alias{testPairedEndBam,BamFile-method}

\title{
    Quickly test if a BAM file has paired end reads
}

\description{
    Quickly test if a BAM file has paired end reads 
}

\usage{
    testPairedEndBam(file, index, ..., fast=TRUE)
    \S4method{testPairedEndBam}{character}(file,index, ..., fast=TRUE)
    \S4method{testPairedEndBam}{BamFile}(file, index, ..., fast=TRUE)
}

\arguments{
    \item{file}{The character(1) file name of the `BAM' ('SAM' for
        \code{asBam}) file to be processed.}
        
    \item{index}{The character(1) name of the index file of the 'BAM' file
        being processed; this is given \emph{without} the '.bai' extension.}
    
    \item{\dots}{Additional arguments, passed to methods.}
    
    \item{fast}{If fast is TRUE(default), then 1000000 reads are used to 
        determine if the BAM file has paired end reads, otherwise all reads of 
        the BAM file are used.}
}

\value{
    A logical vector of length 1 containing TRUE is returned if BAM file 
    contained paired end reads, FALSE otherwise.  
}

\author{
    Martin Morgan <mtmorgan@fhcrc.org>, Sonali Arora <sarora@fhcrc.org>
}

\examples{
    fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
    testPairedEndBam(fl)
}

