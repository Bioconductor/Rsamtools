# Migration of Rsamtools to Rhtslib


## What is the migration about?

Before this migration (i.e. for versions < 1.99.0), Rsamtools contained a
copy of an old (9-10 year old?) version of the samtools and tabix C code.

The latest version of samtools (as of Feb 6, 2019) is 1.9:
http://www.htslib.org/

In the recent years the samtools code has been split into samtools + htslib.
Most of the code that used to be in samtools is now in htslib. The latest
version of htslib is also 1.9.

Bioconductor package Rhtslib contains htslib 1.7, which is a decently recent
version of htslib.

After this migration (i.e. for versions >= 1.99.0), Rsamtools no longer
contains a copy of the samtools and tabix C code but compiles and links
against Rhtslib instead (i.e. against htslib 1.7).


## Current migration status

### C/C++ code

In Rsamtools 1.99.0 all the C/C++ code from the old Rsamtools was migrated to
Rhtslib, except for the code behind `asBcf()` and `razip()`:

- `asBcf()` (implemented in `src/bcffile.c`) is currently disabled but will
  need to be migrated if it turns out that some other package needs it (it
  doesn't seem to be the case though). Note that `asBcf()` isn't used anywhere
  in Rsamtools either (no example in the man page, no unit test, not used in
  the vignette).

- `razip()` (implemented in `zip_compression.c` in the old Rsamtools) was
  removed. This is because recent samtools and htslib no longer support RAZF.
  From the samtools NEWS file:

    RAZF and razip are superseded by BGZF/bgzip and have been removed from
    samtools.

  This happened in samtools 1.0 (15 August, 2014).

### Other remaining problems

- Even though the C code behind `applyPileups()` was migrated (`pileupbam.c`),
  the examples in its man page and most of the tests in
  `inst/unitTests/test_applyPileups.R` are currently failing and have been
  disabled for now.

- Test `test_BcfFile_scanBcfHeader_remote` (located in
  `inst/unitTests/test_BcfFile.R`) passes when run interactively (or via
  `BiocGenerics:::testPackage('Rsamtools')`) but, for some obscure reason,
  NOT in the context of `R CMD check`. It has been disabled for now.


## Status of Bioconductor packages that depend directly on Rsamtools

As of Feb 6, 2019, 159 Bioconductor packages (154 software, 3 data-experiment,
and 2 workflow packages) depend **directly** on Rsamtools (i.e. via their
Depends, Imports, or LinkingTo field). So we've only checked manually a few
of them. This was on a 64-bit Ubuntu 16.04.5 LTS laptop running R 3.6
(2018-12-04 r75758) and Bioconductor 3.9 (current devel).

### All software packages with "LinkingTo: Rsamtools" were tested

  - VariantAnnotation: migrated (changes committed),
    passes `R CMD check`

  - ArrayExpressHTS: migrated (changes not committed yet, Angela Goncalves
    contacted but didn't answer yet),
    passes `R CMD check`

  - BitSeq: migrated (changes committed),
    passes `R CMD check`

  - DiffBind: migrated (changes committed by Gord Brown),
    passes `R CMD check`

  - h5vc: migrated (changes committed, including fixing pre-existing
    build ERROR currently visible on the build report),
    passes `R CMD check`

  - podkat: migrated (changes committed, including fixing
    `inst/examples/example1.vcf.gz` to make it compatible with bcftools 1.7),
    passes `R CMD check`

  - qrqc: migrated (PR created at https://github.com/vsbuffalo/qrqc/pull/6
    per Vince Buffalo's request), passes `R CMD check`

  - QuasR: migrated (PR created at https://github.com/fmicompbio/QuasR/pull/9
    per Michael Stadler's request), passes `R CMD check`

  - seqbias: migrated (patch sent to, and applied by, Daniel Jones),
    passes `R CMD check`

  - TransView: migrated (changes committed),
    passes `R CMD check`

### All software packages that use applyPileups() were tested

  - AllelicImbalance: migrated (contained old BCF file (`ERP000101.bcf`, in
    `extdata/ERP000101_subset/`) that needed to be fixed and reindexed
    because it no longer worked with recent bcftools/Rhtslib)
    passes `R CMD check`

  - biovizBase: passes `R CMD check` _as-is_

  - compEpiTools: passes `R CMD check` _as-is_

  - SICtools: `test_snpDiff()` (from `test_snpDiff.R`) FAILS!

  - VariantTools: passes `R CMD check` _as-is_

### A random sample of a few other software packages were tested

  - AllelicImbalance: migrated (changes not committed yet)
    passes `R CMD check`

  - BaalChIP: passes `R CMD check` _as-is_

  - chimera: passes `R CMD check` _as-is_

  - CoverageView: passes `R CMD check` _as-is_

  - exomeCopy: passes `R CMD check` _as-is_

  - exomePeak: passes `R CMD check` _as-is_

  - GenomicAlignments: passes `R CMD check` _as-is_

  - gmapR: does NOT compile (needs to be migrated)

  - MEDIPS: passes `R CMD check` _as-is_

  - methylPipe: passes `R CMD check` _as-is_

  - rnaSeqMap: passes `R CMD check` _as-is_

  - ssviz: passes `R CMD check` _as-is_

  - systemPipeR: passes `R CMD check` _as-is_

### All workflow packages were tested

  - rnaseqGene: vignette builds _as-is_

  - sequencing: vignette does NOT build (because of AnnotationHub razip'ed
    FASTA file `AH18522`)


## Status of CRAN packages that depend directly on Rsamtools

As of Feb 6, 2019, 8 CRAN packages depend **directly** on Rsamtools (i.e.
via their Depends, Imports, or LinkingTo field).

- BinQuasi: passes `R CMD check` _as-is_

- Brundle: passes `R CMD check` _as-is_

- ExomeDepth: passes `R CMD check` _as-is_

- hoardeR: fails to install because of unrelated pre-existing problem:
  ```
  Error : object ‘importGFF3’ is not exported by 'namespace:GenomicTools'
  ```
  (See https://www.r-project.org/nosvn/R.check/r-devel-linux-x86_64-debian-gcc/hoardeR-00install.html)

- NIPTeR: passes `R CMD check` _as-is_

- PlasmaMutationDetector: passes `R CMD check` _as-is_

- RAPIDR: passes `R CMD check` _as-is_

- spp: passes `R CMD check` _as-is_


## What to do about razip-compressed FASTA files and their index files?

### The problem

RAZF is no longer supported in recent samtools, and old razip-compressed
FASTA files (`.rz`) and their index files (`rz.fai`) are now causing problems.
Here is an example from the sequencing workflow:
```
library(AnnotationHub)
ah <- AnnotationHub()
fa <- ah[["AH18522"]]
library(Rsamtools)
idx <- scanFaIndex(fa)  # still works
long <- idx[width(idx) > 82000]
getSeq(fa, param=long)  # ERROR! ('open' index failed)
```

### What to do about it?

The new way to go is to compress with bgzip and to recreate the index.

#### From the Unix command line

```
/path/to/samtools-1.7/htslib/bgzip myfile.fa       # creates myfile.fa,gz

/path/to/samtools-1.7/samtools faidx myfile.fa,gz  # generates index files
                                                   # myfile.fa.gz.fai
                                                   # and myfile.fa.gz.gzi
```

Note that the compression step is actually optional i.e. one can use
`samtools faidx` directly on the uncompressed FASTA file:

```
gunzip myfile.fa,gz  # uncompress to get myfile.fa back
/path/to/samtools-1.7/samtools faidx myfile.fa     # generates index file
                                                   # myfile.fa.fai only
```

Note that in this case, only the `.fai` index file is generated
(no `.gzi` index file).

#### From R

```
library(Rsamtools)
bgzip("myfile.fa")
indexFa("myfile.fa.bgz")
```

Then create a FaFile object that can be used with `getSeq()` as usual:
```
fa <- FaFile("myfile.fa.bgz")
```


## TODO

Here is an itemized summary of what still needs to happen (some details are
provided above in this document):

- Chase down old razip'ed FASTA files and update them. Places to look at:
  * AnnotationHub
  * ExperimentHub
  * data annotation packages
  * data experiment packages
  * software packages

- Troubleshoot `applyPileups()` unit test failures.

- Troubleshoot `test_BcfFile_scanBcfHeader_remote()` failure (test is
  located in `inst/unitTests/test_BcfFile.R`).

- Migrate `asBcf()`.

- Update `NEWS` file.

