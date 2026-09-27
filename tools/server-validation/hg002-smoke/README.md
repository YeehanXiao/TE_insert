# HG002 ONT smoke test

This directory provides a compact, real-genome smoke test for TEi's experimental
long-read workflow. The preparation script retrieves only a 520,001 bp interval
from an indexed Genome in a Bottle (GIAB) HG002 Oxford Nanopore genomic-DNA BAM;
it does not download the 174.64 GiB source BAM.

The fixed interval is `chr1:74650000-75170000` on GRCh38. It lies completely
inside the GIAB v5.0q structural-variant benchmark interval
`chr1:74378646-75760494` and contains two sequence-resolved insertion records:

| VCF position | Inserted length | Genotype | RepeatMasker annotation |
| --- | ---: | --- | --- |
| `chr1:74727205` | 3,710 bp | `0|1` | `L1HS`, `LINE/L1` |
| `chr1:75091476` | 314 bp | `1|1` | `ALUYA5`, `SINE/Alu` |

`expected_loci.tsv` records the one-based VCF positions and the corresponding
zero-based interbase breakpoints used by TEi. For these simple insertion alleles,
the insertion follows the VCF REF anchor, so `expected_breakpoint_0` equals the
VCF POS value.

## Files

- `prepare_hg002_smoke_test.sh`: downloads the source BAI to a temporary
  directory and uses HTTP range requests through samtools to extract the fixed
  region.
- `truth_te_insertions.fa`: the two v5.0q ALT alleles with the shared REF anchor
  removed. Sequence lengths and SHA-256 digests were checked against the source
  VCF.
- `truth_te_metadata.tsv`: TEi metadata keyed by the first FASTA-header token.
- `expected_loci.tsv`: coordinates, genotypes, repeat annotations, benchmark
  interval, and sequence digests.

## Prepare the regional BAM

The Linux server needs `curl` and samtools built with HTTPS support. The script
downloads a 52.10 MiB BAI into a temporary directory. The resulting regional BAM
is normally only tens of MiB, although its exact size depends on the long reads
overlapping the interval.

```bash
cd tools/server-validation/hg002-smoke
chmod +x prepare_hg002_smoke_test.sh
./prepare_hg002_smoke_test.sh ./data 8
```

The output is:

```text
data/HG002_GRCh38_ONT_UL_UCSC_20200508.chr1_74650000_75170000.bam
data/HG002_GRCh38_ONT_UL_UCSC_20200508.chr1_74650000_75170000.bam.bai
```

The script refuses to overwrite an existing regional BAM. Remove or rename the
previous output explicitly before repeating the preparation.

## Run TEi validation

Build the TEi source tarball first, then run the existing isolated server
validation wrapper from this directory:

```bash
../run_validation.sh \
  /path/to/TEi_1.1.0.9000.tar.gz \
  ./data/HG002_GRCh38_ONT_UL_UCSC_20200508.chr1_74650000_75170000.bam \
  ./truth_te_insertions.fa \
  ont \
  ./validation-output \
  ./truth_te_metadata.tsv \
  8
```

Review calls within 100 bp of either expected breakpoint:

```bash
awk -F '\t' '
  NR == 1 || ($2 == "chr1" &&
    (($3 >= 74727105 && $3 <= 74727305) ||
     ($3 >= 75091376 && $3 <= 75091576)))
' validation-output/insertion_calls.tsv
```

An execution that completes with no calls still tests package installation and
data handling, but it does not satisfy the biological positive-control
expectation. Inspect `evidence.tsv`, `classifications.tsv`, and the BAM around the
two loci before changing thresholds.

## Authoritative sources

HG002 ONT genomic-DNA alignment:

- Dataset README, including library preparation, Guppy 3.2.5, minimap2
  `-a -z 600,200 -x map-ont`, WhatsHap phasing, CC0, and the NIST data-use
  policy:
  <https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/README_ONT-UL_UCSC_HG002.md>
- BAM, 187,513,444,960 bytes:
  <https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/HG002_GRCh38_ONT-UL_UCSC_20200508.phased.bam>
- BAI, 54,635,240 bytes:
  <https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/HG002_GRCh38_ONT-UL_UCSC_20200508.phased.bam.bai>

GIAB HG002 v5.0q GRCh38 benchmark:

- NIST GIAB release notice and usage context:
  <https://www.nist.gov/programs-projects/genome-bottle>
- Structural-variant VCF, 48,040,012 bytes:
  <https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_GRCh38_v5.0q_stvar.vcf.gz>
- VCF index, 1,554,436 bytes:
  <https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_GRCh38_v5.0q_stvar.vcf.gz.tbi>
- Benchmark BED, 230,609 bytes:
  <https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_GRCh38_v5.0q_stvar.benchmark.bed>

## Interpretation limits

- GIAB v5.0q is a general assembly-based small- and structural-variant
  benchmark, not a dedicated mobile-element-insertion benchmark. These two
  records are selected because the VCF itself annotates their inserted alleles
  as `L1HS` and `ALUYA5` with RepeatMasker.
- The bundled FASTA contains the exact truth alleles rather than population
  consensus sequences. This intentionally makes classification easier and must
  not be used to estimate real-world TE-family classification accuracy.
- The source reads and the assembly-based benchmark are not an independent
  discovery/validation pair. This bundle is a reproducible positive-control
  smoke test, not a sensitivity or specificity benchmark.
- Regional BAM extraction retains alignments overlapping the interval and their
  SA tags, but supplementary records whose alignments do not overlap the region
  are not copied. TEi also skips hard-clipped primary alignments because their
  full read sequence is unavailable.
- The 2020 alignment uses an older ONT basecaller and minimap2 workflow. It tests
  compatibility with real genomic long reads, not performance on current ONT
  chemistry or current high-accuracy basecalling.
