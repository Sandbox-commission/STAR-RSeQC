#!/usr/bin/env bash
set -euo pipefail

FASTQ_DIR=""
OUTPUT_DIR=""
GENOME_DIR=""
GTF_FILE=""

usage() {
  cat <<USAGE
Usage: scripts/preflight.sh --fastq-dir DIR --output-dir DIR --genome-dir DIR --gtf FILE
USAGE
}

die() {
  echo "[ERROR] $*" >&2
  exit 1
}

parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --fastq-dir)
        FASTQ_DIR="${2:-}"
        shift 2
        ;;
      --output-dir)
        OUTPUT_DIR="${2:-}"
        shift 2
        ;;
      --genome-dir)
        GENOME_DIR="${2:-}"
        shift 2
        ;;
      --gtf)
        GTF_FILE="${2:-}"
        shift 2
        ;;
      -h|--help)
        usage
        exit 0
        ;;
      *)
        die "Unknown option: $1"
        ;;
    esac
  done

  [[ -n "$FASTQ_DIR" ]] || die "--fastq-dir is required"
  [[ -n "$OUTPUT_DIR" ]] || die "--output-dir is required"
  [[ -n "$GENOME_DIR" ]] || die "--genome-dir is required"
  [[ -n "$GTF_FILE" ]] || die "--gtf is required"
}

check_paths() {
  [[ -d "$FASTQ_DIR" ]] || die "FASTQ directory not found: $FASTQ_DIR"
  [[ -d "$GENOME_DIR" ]] || die "Genome index directory not found: $GENOME_DIR"
  [[ -f "$GTF_FILE" ]] || die "GTF file not found: $GTF_FILE"

  mkdir -p "$OUTPUT_DIR"
  [[ -w "$OUTPUT_DIR" ]] || die "Output directory is not writable: $OUTPUT_DIR"
}

check_fastq_pairs() {
  local count=0
  local missing=0

  shopt -s nullglob
  for r1 in "$FASTQ_DIR"/*_1P.fastq.gz; do
    count=$((count + 1))
    local base
    local r2
    base="${r1%_1P.fastq.gz}"
    r2="${base}_2P.fastq.gz"
    if [[ ! -f "$r2" ]]; then
      echo "[ERROR] Missing R2 pair for: $(basename "$r1")" >&2
      missing=1
    fi
  done
  shopt -u nullglob

  [[ $count -gt 0 ]] || die "No *_1P.fastq.gz files found in: $FASTQ_DIR"
  [[ $missing -eq 0 ]] || die "FASTQ pair validation failed"

  echo "[OK] FASTQ pair check passed ($count samples detected)"
}

check_star_index_hint() {
  # STAR index directories usually contain Genome and SA* files.
  local required=(Genome SA SAindex chrLength.txt chrName.txt chrStart.txt)
  local miss=0
  local f
  for f in "${required[@]}"; do
    if [[ ! -e "$GENOME_DIR/$f" ]]; then
      echo "[WARN] STAR index file missing (or custom index layout): $GENOME_DIR/$f"
      miss=1
    fi
  done
  if [[ $miss -eq 0 ]]; then
    echo "[OK] STAR index layout check passed"
  fi
}

main() {
  parse_args "$@"
  check_paths
  check_fastq_pairs
  check_star_index_hint
  echo "[OK] Preflight completed"
}

main "$@"
