#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat >&2 <<'EOF'
Usage:
  run_validation.sh TEi.tar.gz host.bam te_reference.fa PLATFORM OUTPUT_DIR [te_metadata.tsv|-] [THREADS]

PLATFORM must be ont, hifi, or clr. OUTPUT_DIR must not already exist.
Use - when no TE metadata table is supplied. THREADS defaults to 4.
EOF
  exit 2
}

[[ $# -ge 5 && $# -le 7 ]] || usage

for command in R Rscript minimap2 samtools readlink; do
  command -v "$command" >/dev/null 2>&1 || {
    printf 'ERROR: required command not found: %s\n' "$command" >&2
    exit 127
  }
done

[[ -f "$1" ]] || { printf 'ERROR: TEi source tarball not found: %s\n' "$1" >&2; exit 2; }
[[ -f "$2" ]] || { printf 'ERROR: BAM not found: %s\n' "$2" >&2; exit 2; }
[[ -f "$3" ]] || { printf 'ERROR: TE FASTA not found: %s\n' "$3" >&2; exit 2; }

tarball=$(readlink -f "$1")
bam=$(readlink -f "$2")
te_reference=$(readlink -f "$3")
platform=$4
output_dir=$(readlink -m "$5")
te_metadata=${6:--}
threads=${7:-4}

[[ "$platform" =~ ^(ont|hifi|clr)$ ]] || { printf 'ERROR: PLATFORM must be ont, hifi, or clr.\n' >&2; exit 2; }
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { printf 'ERROR: THREADS must be a positive integer.\n' >&2; exit 2; }
[[ ! -e "$output_dir" ]] || { printf 'ERROR: OUTPUT_DIR already exists: %s\n' "$output_dir" >&2; exit 2; }

if [[ "$te_metadata" != "-" ]]; then
  [[ -f "$te_metadata" ]] || { printf 'ERROR: TE metadata TSV not found: %s\n' "$te_metadata" >&2; exit 2; }
  te_metadata=$(readlink -f "$te_metadata")
fi

samtools quickcheck -v "$bam" || { printf 'ERROR: samtools quickcheck failed for %s\n' "$bam" >&2; exit 2; }

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$output_dir/R-library"
export R_LIBS_USER="$output_dir/R-library"
trap 'status=$?; printf "ERROR: validation stopped at line %s with exit status %s.\n" "$LINENO" "$status" >&2; exit "$status"' ERR

{
  printf 'R: '
  R --version | sed -n '1p'
  printf 'minimap2: '
  minimap2 --version
  printf 'samtools: '
  samtools --version | sed -n '1p'
} > "$output_dir/tool_versions.txt"

printf '[1/2] Installing TEi into %s\n' "$R_LIBS_USER"
R CMD INSTALL --preclean --clean --library="$R_LIBS_USER" "$tarball" \
  2>&1 | tee "$output_dir/package_install.log"

printf '[2/2] Running long-read validation\n'
Rscript --vanilla "$script_dir/validate_long_read.R" \
  "$bam" "$te_reference" "$platform" "$output_dir" "$te_metadata" "$threads" \
  2>&1 | tee "$output_dir/validation.log"

printf 'Validation completed: %s\n' "$output_dir"
