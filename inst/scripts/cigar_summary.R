## Simple functions for reading and summarizing CIGARs from BAM files

## Creates a table of CIGAR values from a BAM file
##
## Input:
##   file: The file name of the BAM file to be parsed.
##
## Output:
##   A DataFrame with two columns:
##     cigar (CompressedRleList) and count (integer)
readCigarTable <- function(file)
{
    fileData <- scanBam(file, param = ScanBamParam(what=c("cigar")))
    if (length(fileData) != 1)
        stop("scanBam generated more than 1 list element")
    cigar <-
      factor(fileData[[1]][["cigar"]]@.cigar[!is.na(fileData[[1]][["cigar"]]@.cigar)])
    basicTable <- table(cigar)
    tableOrder <- order(basicTable, decreasing = TRUE)
    cigar <- factor(cigar, levels = levels(cigar)[tableOrder])
    basicTable <- basicTable[tableOrder]
    cigarValues <-
      CharacterList(lapply(strsplit(levels(cigar), "[0-9]+"), "[", -1))
    cigarLengths <- IntegerList(strsplit(levels(cigar), "[A-Za-z]+"))
    DataFrame(cigar =
              IRanges:::newCompressedList("CompressedRleList",
                   Rle(unlist(cigarValues, use.names = FALSE),
                       unlist(cigarLengths, use.names = FALSE)),
                   cumsum(unlist(lapply(cigarLengths, sum)))),
              count = as.integer(basicTable))
}

## Summarize a CIGAR table in the form of a DataFrame
##
## Input:
##   x: a DataFrame with cigar and count columns
##
## Output:
##   A list with two elements:
##     AlignedCharacters (integer) and Indels (matrix)
summarizeCigarTable <- function(x) {
    alignedCharacters <-
      table(rep.int(elementLengths(x[["cigar"]]), x[["count"]]),
            rep.int(viewSums(Views(unlist(x[["cigar"]] == "M"),
                                   as(x[["cigar"]]@partitioning, "IRanges"))) ==
                    elementLengths(x[["cigar"]]), x[["count"]]))
    tabledAlignedCharacters <- as(rev(colSums(alignedCharacters)), "integer")
    names(tabledAlignedCharacters) <-
      unname(c("TRUE" = "AllAligned",
               "FALSE" = "SomeNonAligned")[names(tabledAlignedCharacters)])

    indelHits <-
      rbind(data.frame(subject =
                       subjectHits(findOverlaps(IRanges(unlist(x[["cigar"]] == "D")),
                                                x[["cigar"]]@partitioning)),
                       type = factor("D", levels = c("D", "I"))),
            data.frame(subject =
                       subjectHits(findOverlaps(IRanges(unlist(x[["cigar"]] == "I")),
                                                x[["cigar"]]@partitioning)),
                       type = factor("I", levels = c("D", "I"))))
    tabledIndelHits <- table(indelHits[,1], indelHits[,2])
    tabledIndelHits <-
      tabledIndelHits[rep.int(seq_len(nrow(tabledIndelHits)),
                              x[["count"]][as.integer(rownames(tabledIndelHits))]),
                      , drop = FALSE]
    tabledIndelHits <-
      as(table(tabledIndelHits[,"D"], tabledIndelHits[,"I"]), "matrix")
    rownames(tabledIndelHits) <- paste("D", rownames(tabledIndelHits), sep = "")
    colnames(tabledIndelHits) <- paste("I", colnames(tabledIndelHits), sep = "")

    list("AlignedCharacters" = tabledAlignedCharacters,
         "Indels" = tabledIndelHits)
}


## Examples
if (FALSE) {

suppressMessages(library(Rsamtools))
dataDir <- "/home/biocdev/data_store/1000genomes/extdata"
file1 <- file.path(dataDir, "NA19239.SLX.maq.SRP000033.2009_09.subset.bam")
file2 <- file.path(dataDir, "NA19240.chrom6.SLX.maq.SRP000032.2009_07.subset.bam")
file3 <- file.path(dataDir, "NA19240.chrom1.454.ssaha2.SRP000032.2009_10.bam")

cigarTable1 <- readCigarTable(file1) ## Takes approx  1 minute
cigarTable2 <- readCigarTable(file2) ## Takes approx  1 minute
cigarTable3 <- readCigarTable(file3) ## Takes approx 30 minutes

summarizeCigarTable(cigarTable1)
summarizeCigarTable(cigarTable2)
summarizeCigarTable(cigarTable3)

}
