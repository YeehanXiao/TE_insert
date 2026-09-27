#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat >&2 <<'EOF'
Usage: prepare_hg002_smoke_test.sh [OUTPUT_DIR] [THREADS]

OUTPUT_DIR defaults to ./data. THREADS defaults to 4.
The output BAM must not already exist.
EOF
  exit 2
}

[[ $# -le 2 ]] || usage

output_dir=${1:-./data}
threads=${2:-4}
region='chr1:74650000-75170000'
bam_name='HG002_GRCh38_ONT_UL_UCSC_20200508.chr1_74650000_75170000.bam'
source_bam='https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/HG002_GRCh38_ONT-UL_UCSC_20200508.phased.bam'
source_bai="${source_bam}.bai"

[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { printf 'ERROR: THREADS must be a positive integer.\n' >&2; exit 2; }
for command in awk curl mktemp samtools; do
  command -v "$command" >/dev/null 2>&1 || {
    printf 'ERROR: required command not found: %s\n' "$command" >&2
    exit 127
  }
done

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
awk '
  /^>/ {
    if (id != "") length_by_id[id] = length(sequence)
    id = substr($1, 2)
    sequence = ""
    next
  }
  $0 !~ /^[ACGTN]+$/ { invalid = 1 }
  { sequence = sequence $0 }
  END {
    if (id != "") length_by_id[id] = length(sequence)
    exit(invalid ||
      length_by_id["L1HS_HG002_GRCh38_chr1_74727205_GIAB_v5.0q"] != 3710 ||
      length_by_id["ALUYA5_HG002_GRCh38_chr1_75091476_GIAB_v5.0q"] != 314)
  }
' "$script_dir/truth_te_insertions.fa" || {
  printf 'ERROR: bundled truth TE FASTA failed validation.\n' >&2
  exit 2
}

mkdir -p "$output_dir"
output_bam="$output_dir/$bam_name"
[[ ! -e "$output_bam" && ! -e "$output_bam.bai" ]] || {
  printf 'ERROR: output already exists: %s\n' "$output_bam" >&2
  exit 2
}

temp_dir=$(mktemp -d "${TMPDIR:-/tmp}/tei-hg002-smoke.XXXXXX")
trap 'rm -rf -- "$temp_dir"' EXIT

printf '[1/3] Downloading the 52.1 MiB source BAM index\n'
curl -L --fail --retry 5 --retry-all-errors \
  -o "$temp_dir/source.bam.bai" "$source_bai"
[[ $(wc -c < "$temp_dir/source.bam.bai") -eq 54635240 ]] || {
  printf 'ERROR: downloaded source BAM index has an unexpected size.\n' >&2
  exit 2
}

printf '[2/3] Extracting %s from the remote 174.64 GiB BAM\n' "$region"
samtools view -@ "$threads" -X -bh \
  -o "$temp_dir/$bam_name" \
  "$source_bam" "$temp_dir/source.bam.bai" "$region"

printf '[3/3] Validating and indexing the regional BAM\n'
samtools quickcheck -v "$temp_dir/$bam_name"
samtools view -H "$temp_dir/$bam_name" | awk -F '\t' '
  $1 == "@SQ" && $2 == "SN:chr1" && $3 == "LN:248956422" { found = 1 }
  END { exit(!found) }
' || { printf 'ERROR: regional BAM is not aligned to the expected GRCh38 chr1.\n' >&2; exit 2; }
samtools index -@ "$threads" "$temp_dir/$bam_name"
read_count=$(samtools view -c "$temp_dir/$bam_name")
[[ "$read_count" -gt 0 ]] || { printf 'ERROR: regional BAM contains no alignments.\n' >&2; exit 2; }
mv "$temp_dir/$bam_name" "$temp_dir/$bam_name.bai" "$output_dir/"

printf 'Prepared: %s\nRegion: %s\nAlignments: %s\n' "$output_bam" "$region" "$read_count"
