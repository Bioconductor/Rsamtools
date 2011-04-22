.MpileupParam <-
    function(## limit how reads counted
             ## ceilingMapQuality=60L,
             flag=scanBamFlag(),
             minBaseQuality=13L,
             minMapQuality=0L,
             minDepth=0L,
             maxDepth=250L,
             maxIndelSampleDepth=250L,
             minIndelReadCount=1L,
             minIndelReadFraction=0.0002,
             ## coefficients
             gapOpenErrPr=40L,
             gapExtendErrPr=20L,
             homopolymerErrCoef=100,
             useBAQ=TRUE,
             useExtendedBAQ=FALSE,
             ## output
             genotypeLikelihoods=FALSE,
             perSampleReadDepth=FALSE,
             perSampleStrandBiasPr=FALSE,
             noIndels=FALSE,
             ## select reads
             usePlatformsForIndels=character(),
             excludeReadGroups=character(),
             useAnomalousReadPairs=FALSE,
             which=GRanges())
{
    as.list(environment())
}

setGeneric(".asSpace", function(x) standardGeneric(".asSpace"))

setMethod(.asSpace, "RangesList", function(x) {
    list(as.character(space(x)), .uunlist(start(x)), .uunlist(end(x)))
})

setMethod(.asSpace, "GRanges", function(x) {
    list(as.character(seqnames(x)), start(x), end(x))
})

.mpileupBam <-
    function(files, callback=identity, ..., param)
{
    tryCatch({
        files <- lapply(files, .extptr)
        param <- as.list(param)
        space <-
            if (0L != length(param[["which"]])) .asSpace(param[["which"]])
            else NULL
        .Call(.mpileup_bam, files, space, param, callback)
    }, error=function(err) {
        stop("mpileupBam: ", conditionMessage(err), call.=FALSE)
    })
}

## setGeneric("mpileupBam",
##     function(file, index=file, ..., param=MpileupParam())
##         standardGeneric("mpileupBam"),
##     signature="file")


## Usage:   samtools mpileup [options] in1.bam [in2.bam [...]]

## Options: -f FILE     reference sequence file [null]
##          -r STR      region in which pileup is generated [null]
##          -l FILE     list of positions (format: chr pos) [null]
##          -b FILE     list of input BAM files [null]
##          -M INT      cap mapping quality at INT [60]
##          -Q INT      min base quality [13]
##          -q INT      filter out alignment with MQ smaller than INT [0]
##          -d INT      max per-BAM depth [250]
##          -L INT      max per-sample depth for INDEL calling [250]
##          -P STR      comma separated list of platforms for indels [all]
##          -o INT      Phred-scaled gap open sequencing error probability [40]
##          -e INT      Phred-scaled gap extension seq error probability [20]
##          -h INT      coefficient for homopolyer errors [100]
##          -m INT      minimum gapped reads for indel candidates [1]
##          -F FLOAT    minimum fraction of gapped reads for candidates [0.002]
##          -G FILE     exclude read groups listed in FILE [null]
##          -A          use anomalous read pairs in SNP/INDEL calling
##          -g          generate BCF output
##          -u          do not compress BCF output
##          -E          extended BAQ for higher sensitivity but lower specificity
##          -B          disable BAQ computation
##          -D          output per-sample DP
##          -S          output per-sample SP (strand bias P-value, slow)
##          -I          do not perform indel calling

