### test_readGAlignmentPairsFromBam.R

# Flag bits summary
# -----------------
#   0x1: template having multiple segments in sequencing
#   0x2: each segment properly aligned according to the aligner
#   0x4: segment unmapped
#   0x8: next segment in the template unmapped
#  0x10: SEQ being reverse complemented
#  0x20: SEQ of the next segment in the template being reversed
#  0x40: the first segment in the template
#  0x80: the last segment in the template
# 0x100: secondary alignment
# 0x200: not passing quality controls
# 0x400: PCR or optical duplicate

make_samline <- function(QNAME, mapped=0, RNAME="*", POS=0, strand="+",
                         primary=1, flagbits=0,
                         RNEXT="*", PNEXT=0, proper=0)
{
    FLAG <- 0x4 * (!mapped) + 0x10 * (strand == "-") +
            0x100 * (!primary) + 0x2 * proper + flagbits
    if (mapped) {
        CIGAR <- "14M"
    } else {
        CIGAR <- "*"
    }
    paste(QNAME, FLAG, RNAME, POS, "255", CIGAR,
          RNEXT, PNEXT, "0", "*", "*", sep="\t")
}

make_mapped_pair <- function(QNAME,
                             RNAME1, POS1, strand1,
                             RNAME2, POS2, strand2,
                             primary, proper,
                             pair_id=NA)
{
    flag0 <- 0x1 + 0x2 * proper + 0x100 * (!primary)
    FLAG1 <- flag0 + 0x40 + 0x10 * (strand1 == "-") +
                            0x20 * (strand2 == "-")
    FLAG2 <- flag0 + 0x80 + 0x10 * (strand2 == "-") +
                            0x20 * (strand1 == "-")
    line1 <- paste(QNAME, FLAG1, RNAME1, POS1, "255", "18M",
                   RNAME2, POS2, "0", "*", "*", sep="\t")
    line2 <- paste(QNAME, FLAG2, RNAME2, POS2, "255", "18M",
                   RNAME1, POS1, "0", "*", "*", sep="\t")
    if (!identical(pair_id, NA)) {
        pair_id_tag <- paste0("pi:Z:", pair_id)
        line1 <- paste(line1, pair_id_tag, sep="\t")
        line2 <- paste(line2, pair_id_tag, sep="\t")
    }
    c(line1, line2)
}

make_mapped_pairs <- function(mapped_pair_table)
{
    if (!is.data.frame(mapped_pair_table))
        mapped_pair_table <- read.table(mapped_pair_table, header=TRUE,
                                        stringsAsFactors=FALSE)
    row_groups <- unname(split(seq_len(nrow(mapped_pair_table)),
                               mapped_pair_table$QNAME))

    unlist(lapply(row_groups,
        function(row_group) {
            lines <- unlist(lapply(row_group, function(i)
                do.call("make_mapped_pair", mapped_pair_table[i, ])))
            if (length(row_group) == 3L) {
                lines[c(2L, 4L, 6L)] <- lines[c(6L, 2L, 4L)]
            } else {
                lines[c(FALSE, TRUE)] <- rev(lines[c(FALSE, TRUE)])
            }
            lines
        }))
}

make_toy_bamfile <- function(mapped_pair_table, filename)
{
    lines0 <- c("@HD\tVN:1.3",
                "@SQ\tSN:chr1\tLN:2450",
                "@SQ\tSN:chr2\tLN:1882",
                "@SQ\tSN:chrX\tLN:999")
    ## Single end reads
    lines1 <- c(
        ## s001: 1 primary alignment
        make_samline("s001", mapped=1, RNAME="chr1", POS=10, strand="+",
                     primary=1),
        ## s002: 1 primary alignment + 3 secondary alignments
        make_samline("s002", mapped=1, RNAME="chr1", POS=20, strand="+",
                     primary=1),
        make_samline("s002", mapped=1, RNAME="chr1", POS=21, strand="+",
                     primary=0),
        make_samline("s002", mapped=1, RNAME="chr1", POS=22, strand="+",
                     primary=0),
        make_samline("s002", mapped=1, RNAME="chr1", POS=20, strand="+",
                     primary=0),
        ## s003: unmapped
        make_samline("s003")
    )
    ## Paired end reads
    lines2 <- c(
        ## Mapped pairs
        make_mapped_pairs(mapped_pair_table),
        ## p991: 1 pair with a missing mate (can happen if file was subsetted
        ## with e.g. filterBam)
        make_mapped_pair("p991",
                         "chr2", 150, "+", "chr2", 199, "-",
                         primary=1, proper=1)[2L],
        ## p992: 1 pair with 1st mate unmapped and 2nd mate mapped
        make_samline("p992", flagbits=0x1 + 0x40,
                             RNEXT="chr2", PNEXT=150),
        make_samline("p992", mapped=1, RNAME="chr2", POS=150, strand="+",
                             primary=1, flagbits=0x1 + 0x8 + 0x80),
        ## p993: 1 pair with both mates unmapped
        make_samline("p993", flagbits=0x1 + 0x8 + 0x40),
        make_samline("p993", flagbits=0x1 + 0x8 + 0x80)
    )
    ## Reads with multiple segments
    lines3 <- c(
        ## m001: 3 segments in the template (index of each segment in template
        ## is known)
        make_samline("m001", mapped=1,, RNAME="chrX", POS=10, strand="+",
                     flagbits=0x1 + 0x40,
                     RNEXT="chrX", PNEXT=20, proper=1),
        make_samline("m001", mapped=1, RNAME="chrX", POS=20, strand="+",
                     flagbits=0x1 + 0x40 + 0x80,
                     RNEXT="chrX", PNEXT=30, proper=1),
        make_samline("m001", mapped=1, RNAME="chrX", POS=30, strand="+",
                     flagbits=0x1 + 0x80,
                     RNEXT="chrX", PNEXT=10, proper=1),
        ## m002: 3 segments in the template (index of each segment in template
        ## was lost)
        make_samline("m002", mapped=1, RNAME="chrX", POS=10, strand="+",
                     flagbits=0x1,
                     RNEXT="chrX", PNEXT=20, proper=1),
        make_samline("m002", mapped=1, RNAME="chrX", POS=20, strand="+",
                     flagbits=0x1,
                     RNEXT="chrX", PNEXT=30, proper=1),
        make_samline("m002", mapped=1, RNAME="chrX", POS=30, strand="+",
                     flagbits=0x1,
                     RNEXT="chrX", PNEXT=10, proper=1)
    )
    samfile <- paste0(filename, ".sam")
    cat(c(lines0, lines1, lines2, lines3), file=samfile, sep="\n")
    bamfile <- asBam(samfile, filename, overwrite=TRUE)
    ## Should never happen.
    if (bamfile != paste0(filename, ".bam"))
        stop("asBam() returned an unexpected path")
    bamfile
}

### 1 line per mapped pair. Each line will generate 2 lines/records in the
### SAM file.
mapped_pair_table <- textConnection("
QNAME   RNAME1  POS1    strand1 RNAME2  POS2    strand2 primary proper  pair_id
# p001: 1 primary proper pair
p001    chr2    10      +       chr2    110     -       1       1       p001
# p002: 1 primary non proper pair
p002    chr2    20      +       chr2    120     -       1       0       p002
# p003: 2 proper pairs: 1 primary + 1 secondary
p003    chr2    30      +       chr2    130     -       1       1       p003a
p003    chr2    31      +       chr2    131     -       0       1       p003b
# p004: 2 non proper pairs: 1 primary + 1 secondary
p004    chr2    40      +       chr2    140     -       1       0       p004a
p004    chr2    41      +       chr2    141     -       0       0       p004b
# p005: 2 primary pairs (some aligners seem to produce that, even though they
# probably shouldn't)
p005    chr2    50      +       chr2    150     -       1       1       p005a
p005    chr2    51      +       chr2    151     -       1       1       p005b
# p006: 3 pairs: 1 primary proper + 1 secondary proper + 1 secondary non
# proper
p006    chr2    60      +       chr2    160     -       1       1       p006a
p006    chr2    61      +       chr2    161     -       0       1       p006b
p006    chr2    62      +       chr2    60      -       0       0       p006c
# p007: 2 pairs mapped to the same position: 1 primary proper + 1 secondary
# proper
p007    chr2    70      +       chr2    170     -       1       1       p007a
p007    chr2    70      +       chr2    170     -       0       1       p007b
# p008: 3 pairs mapped to the same position: 1 primary proper + 1 secondary
# proper + 1 secondary non proper
p008    chr2    80      +       chr2    180     -       1       1       p008a
p008    chr2    80      +       chr2    180     -       0       1       p008b
p008    chr2    80      +       chr2    180     -       0       0       p008c
# p009: 3 pairs mapped to the same position: 1 primary proper + 2 secondary
# proper. The secondary pairs can NOT be disambiguated.
p009    chr2    90      +       chr2    190     -       1       1       p009a
p009    chr2    90      +       chr2    190     -       0       1       p009b
p009    chr2    90      +       chr2    190     -       0       1       p009c
")

toy_bamfile <- make_toy_bamfile(mapped_pair_table, tempfile())

test_readGAlignmentPairsFromBam <- function()
{
    param <- ScanBamParam(tag="pi")
    galp <- suppressWarnings(
        readGAlignmentPairsFromBam(toy_bamfile, use.names=TRUE, param=param)
    )

    ## Check the dumped alignments
    dumped_gal <- getDumpedAlignments()
    checkTrue(validObject(dumped_gal, complete=TRUE))
    pi_target <- rep(c("p009b", "p009c"), each=2)
    checkIdentical(pi_target, sort(mcols(dumped_gal)$pi))

    ## Check 'galp'
    checkTrue(validObject(galp, complete=TRUE))
    pi_target <- c("p001", "p002", "p003a", "p003b", "p004a", "p004b",
                   "p005a", "p005b", "p006a", "p006b", "p006c",
                   "p007a", "p007b", "p008a", "p008b", "p008c", "p009a")
    checkIdentical(pi_target, mcols(first(galp))$pi)
    checkIdentical(pi_target, mcols(last(galp))$pi)
}

test_readGAlignmentPairsFromBam_which <- function()
{
    ## 4 non-overlapping regions of interest: first two regions only overlap
    ## with first p001 mate and last two regions only with last p001 mate.
    my_ROI <- GRanges("chr2", IRanges(c(10, 15, 110, 115), width=1))
    my_ROI_labels <- c("chr2:10-10", "chr2:15-15",
                       "chr2:110-110", "chr2:115-115")
    param <- ScanBamParam(tag="pi", which=my_ROI[c(1, 4)])
    target1 <- readGAlignmentPairsFromBam(toy_bamfile, use.names=TRUE,
                                          param=param, with.which_label=TRUE)
    checkTrue(validObject(target1, complete=TRUE))
    checkIdentical(1L, length(target1))
    checkIdentical(Rle(factor(my_ROI_labels[1], levels=my_ROI_labels[c(1, 4)])),
                   mcols(first(target1))$which_label)
    checkIdentical(Rle(factor(my_ROI_labels[4], levels=my_ROI_labels[c(1, 4)])),
                   mcols(last(target1))$which_label)
    mcols(target1@first)$which_label <- mcols(target1@last)$which_label <- NULL

    ## Checking all possible combinations of ranges in 'which'.
    check_my_ROI_subsets <- function(subset)
    {
        check_my_ROI_subset <- function(i)
        {
            #print(i)
            param <- ScanBamParam(tag="pi", which=my_ROI[i])
            current <- suppressWarnings(
                readGAlignmentPairsFromBam(toy_bamfile, use.names=TRUE,
                                           param=param)
            )
            if (sum(i <= 2L) == 1L && sum(i >= 3L) == 1L) {
                checkIdentical(target1, current)
            } else {
                checkIdentical(0L, length(current))
            }
        }
        check_my_ROI_subset(subset)
        if (length(subset) >= 2L) {
            check_my_ROI_subset(rev(subset))
            if (length(subset) >= 4L)
                check_my_ROI_subset(c(4L, 1:3))
        }
        TRUE
    }
    for (m in 1:length(my_ROI))
        combn(length(my_ROI), m, FUN=check_my_ROI_subsets)
}

