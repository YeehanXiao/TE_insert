# TEi long-read server validation

This bundle validates a built TEi source tarball against a real host-aligned long-read BAM on Linux. It installs TEi into a private R library inside the output directory and does not install or modify system dependencies.

## Required inputs

- A TEi source tarball produced by `R CMD build`.
- A host-aligned BAM containing read sequences, CIGAR strings, and SA tags for split-alignment evidence.
- A TE FASTA. The first token of each FASTA header is used as the TE identifier.
- One platform name: `ont`, `hifi`, or `clr`.
- An optional tab-separated metadata file with `te_id`, `te_family`, and `te_subfamily` columns.

The server must already provide R 4.2 or newer with the TEi package dependencies, minimap2, samtools, and a working C++ compiler. The script stops with a nonzero exit status when a required command, input, R dependency, package build step, or analysis step fails.

## Run

Without TE metadata:

```bash
chmod +x run_validation.sh
./run_validation.sh \
  ./TEi_1.1.0.9000.tar.gz \
  /data/sample.host.bam \
  /data/te_reference.fa \
  ont \
  ./validation_ont \
  - \
  16
```

With TE metadata:

```bash
./run_validation.sh \
  ./TEi_1.1.0.9000.tar.gz \
  /data/sample.host.bam \
  /data/te_reference.fa \
  hifi \
  ./validation_hifi \
  /data/te_metadata.tsv \
  16
```

Use a new output directory for every run. The default thread count is 4 when the final argument is omitted.

## HG002 ONT positive-control smoke test

The bundled `hg002-smoke` kit extracts only `chr1:74650000-75170000` from the
official GIAB HG002 ONT ultra-long genomic-DNA BAM. The 520 kb interval contains
sequence-resolved L1HS and AluYa5 truth insertions. It avoids downloading the
complete 174.64 GiB source BAM.

On the Linux server, run:

```bash
cd hg002-smoke
chmod +x prepare_hg002_smoke_test.sh ../run_validation.sh
./prepare_hg002_smoke_test.sh ./data 8
../run_validation.sh \
  ../TEi_1.1.0.9000.tar.gz \
  ./data/HG002_GRCh38_ONT_UL_UCSC_20200508.chr1_74650000_75170000.bam \
  ./truth_te_insertions.fa \
  ont \
  ./validation-output \
  ./truth_te_metadata.tsv \
  8
```

See [`hg002-smoke/README.md`](hg002-smoke/README.md) for the expected loci,
source records, call-review command, and interpretation limits.

## Fixed validation settings

The validation uses host mapping quality 20, minimum candidate length 100 bp, maximum host reference gap 50 bp, minimum TE-aligned length 80 bp, minimum TE query fraction 0.5, minimum TE identity 0.7, clustering window 50 bp, and minimum support of two unique reads. Ambiguous TE assignments are retained in the classification output but excluded from insertion calls.

TE classification uses the corresponding minimap2 preset: `map-ont`, `map-hifi`, or `map-pb` for ONT, HiFi, or CLR data.

## Outputs

- `validation_summary.tsv`: package version and record counts.
- `evidence.tsv`: extracted CIGAR insertion, terminal soft-clip, and split-gap evidence.
- `classifications.tsv`: best TE assignment and ambiguity status for each evidence record.
- `te_hits.tsv`: all parsed minimap2 hits.
- `insertion_calls.tsv`: clustered insertion calls and support filters.
- `parameters.tsv`: exact analysis parameters used for the run.
- `long_read_result.rds`: complete R result object.
- `limitations.txt`: current experimental limitations.
- `sessionInfo.txt` and `tool_versions.txt`: reproducibility information.
- `package_install.log` and `validation.log`: complete logs.
- `R-library/`: isolated library containing the tested TEi installation.

A successful run with zero insertion calls is still a valid execution result. Review `validation_summary.tsv`, evidence counts, alignment quality, TE reference content, and thresholds before drawing a biological conclusion.
