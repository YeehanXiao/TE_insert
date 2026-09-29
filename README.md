# TEi

[![R-CMD-check](https://github.com/YeehanXiao/TE_insert/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/YeehanXiao/TE_insert/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE.md)

TEi detects candidate non-reference transposable element insertions from
genome-aligned sequencing reads. The maintained package contains two explicit
workflows:

- a short-read workflow based on terminal soft clips and TE-reference
  realignment;
- an experimental long-read workflow for Oxford Nanopore and PacBio BAM files.

The long-read workflow reports evidence and ambiguity; it does not claim
genotyping, local assembly, mosaic calibration, or complete insertion
reconstruction.

## Installation

Install the required Bioconductor packages and then TEi:

```r
install.packages("BiocManager")
BiocManager::install(c(
  "Biostrings", "GenomicRanges", "IRanges", "Rhtslib",
  "Rsamtools", "S4Vectors"
))
install.packages("remotes")
remotes::install_github("YeehanXiao/TE_insert")
```

The short-read TE-alignment step additionally uses `Rbowtie2`:

```r
BiocManager::install("Rbowtie2")
```

Long-read TE classification requires
[minimap2](https://github.com/lh3/minimap2) on `PATH`.

## Short-read workflow

```r
library(TEi)

buildIndex("te_consensus.fa", "te_index")

extractSoftClip(
  file = "host.qname_sorted.bam",
  outfq = "soft_clips.fastq",
  mapq = 20,
  length = 30,
  tsd = 10
)

te_bam <- alignment(
  reference = "te_index",
  fastq = "soft_clips.fastq",
  bamOutput = "soft_clips_to_te.bam",
  threads = 8
)

insertLocation(te_bam, "te_breakpoints.bed", ratio = 0.5)
calls <- processInsertion("te_breakpoints.bed", max.gapwidth = 10)
```

The TE-aligned BAM must be query-name sorted. `alignment()` performs that step
automatically.

A complete WDL 1.0 implementation of the same three-stage workflow is installed
with the package:

```r
system.file(
  "workflows", "short-read", "TEi-short-read.wdl",
  package = "TEi", mustWork = TRUE
)
system.file(
  "workflows", "short-read", "inputs.example.json",
  package = "TEi", mustWork = TRUE
)
```

Copy the example inputs file, replace its two input paths, and run it with
Cromwell or another WDL 1.0 runner. The workflow directly calls the package
functions and does not require separate wrapper scripts.

## Experimental long-read workflow

Input must be genomic DNA reads aligned to a host reference in BAM format.
Recommended host-alignment presets are `map-ont` for standard ONT reads,
`lr:hq` for accurate ONT reads, `map-hifi` for PacBio HiFi, and `map-pb` for
legacy PacBio CLR data. Preserve soft clipping and SA tags in the host BAM.

```r
evidence <- extractLongReadEvidence(
  alignment = "long_reads.host.bam",
  minMapQ = 20,
  minLength = 100
)

classification <- classifyLongReadEvidence(
  evidence = evidence,
  teReference = "te_consensus.fa",
  teMetadata = "te_metadata.tsv",
  platform = "ont",
  threads = 8
)

calls <- summarizeLongReadInsertions(
  evidence,
  classification,
  clusterWindow = 50,
  minSupport = 2
)
```

The convenience wrapper runs all three stages:

```r
result <- detectLongReadInsertions(
  alignment = "long_reads.host.bam",
  teReference = "te_consensus.fa",
  platform = "hifi",
  threads = 8
)
```

Long-read evidence currently includes large CIGAR insertions, terminal soft
clips, and colinear same-strand gaps reconstructed from SA tags. Candidate
segments are mapped to a user-supplied TE FASTA with minimap2. Near-tied hits
remain `ambiguous` rather than being forced into one TE family.

The experimental long-read workflow extends TEi's original short-read
evidence--classification--aggregation design to ONT and PacBio alignments. Its
implementation was newly written for this repository using Rsamtools and
minimap2; it does not reuse source code from third-party TE-insertion callers.

TE reference sequences are not bundled. Users are responsible for selecting a
licensed and biologically appropriate reference library.

## Output and validation status

Short-read calls are returned as `GRanges`. Long-read functions return explicit
evidence, classification, and call tables with zero-based interbase breakpoint
columns marked by the `_0` suffix.

The package includes synthetic regression tests for native soft-clip handling,
TE multimapping, long-read CIGAR insertions, SA gaps, filtering, ambiguity, and
support aggregation. Real ONT/HiFi sensitivity and precision still require
dataset-specific benchmarking before biological claims are made.

A reproducible positive-control smoke test is provided in
[`tools/server-validation/hg002-smoke`](tools/server-validation/hg002-smoke).
It extracts a 520 kb window from the official GIAB HG002 ONT ultra-long
genomic-DNA alignment and checks two sequence-resolved truth insertions: L1HS
and AluYa5. This is an execution and locus-recovery test, not a sensitivity or
specificity benchmark.

## Contributors and attribution

TEi was originally co-developed and co-authored by
[Yihan Xiao (`@YeehanXiao`)](https://github.com/YeehanXiao) and
[Tao Chen (`@tchen-tt`)](https://github.com/tchen-tt).

- **Yihan Xiao:** original co-developer and package co-author; current
  maintainer and developer of the package modernization, short-read WDL
  workflow, test infrastructure, and experimental long-read extension.
- **Tao Chen:** original co-developer and package co-author.

This repository retains the complete commit history from the original
[`tchen-tt/TEi`](https://github.com/tchen-tt/TEi) project and its MIT license
notice. Current maintenance focuses on reproducible packaging, tests, and
explicit short-read and experimental long-read interfaces.

## References

- Li H. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*
  (2018).
- Tang Z et al. Human transposon insertion profiling. *PNAS* (2017).
- Chu C et al. Comprehensive identification of transposable element insertions
  using multiple sequencing technologies. *Nature Communications* (2021).

## License

MIT. See [LICENSE.md](LICENSE.md).
