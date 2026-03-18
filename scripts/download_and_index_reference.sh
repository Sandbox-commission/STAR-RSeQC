#!/usr/bin/env bash
set -euo pipefail

# Download Ensembl human reference (FASTA + GTF) and build a STAR index.
# Defaults target GRCh38 release 113 to match this pipeline's examples.

RELEASE="113"
ASSEMBLY="GRCh38"
SPECIES="homo_sapiens"
READ_LENGTH="101"
THREADS="$(nproc 2>/dev/null || echo 8)"
OUT_DIR="refs"
STAR_BIN="STAR"

usage() {
  cat <<USAGE
Usage: scripts/download_and_index_reference.sh [options]

Options:
  --release <N>        Ensembl release (default: 113)
  --assembly <NAME>    Assembly name (default: GRCh38)
  --species <NAME>     Species (default: homo_sapiens)
  --read-length <N>    Read length for sjdbOverhang+1 (default: 101)
  --threads <N>        STAR indexing threads (default: nproc)
  --out-dir <DIR>      Output base directory (default: refs)
  --star-bin <PATH>    STAR binary/path (default: STAR)
  -h, --help           Show this help

Output layout (under --out-dir):
  fasta/genome.fa.gz
  gtf/annotation.gtf.gz
  annotation.gtf
  star_index/
USAGE
}

die() {
  echo "[ERROR] $*" >&2
  exit 1
}

have() {
  command -v "$1" >/dev/null 2>&1
}

parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --release)
        RELEASE="${2:-}"
        shift 2
        ;;
      --assembly)
        ASSEMBLY="${2:-}"
        shift 2
        ;;
      --species)
        SPECIES="${2:-}"
        shift 2
        ;;
      --read-length)
        READ_LENGTH="${2:-}"
        shift 2
        ;;
      --threads)
        THREADS="${2:-}"
        shift 2
        ;;
      --out-dir)
        OUT_DIR="${2:-}"
        shift 2
        ;;
      --star-bin)
        STAR_BIN="${2:-}"
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
}

download_file() {
  local url="$1"
  local dest="$2"

  if have curl; then
    curl -fL --retry 3 --retry-delay 2 -o "$dest" "$url"
  elif have wget; then
    wget -O "$dest" "$url"
  else
    die "Neither curl nor wget is installed. Install one and retry."
  fi
}

main() {
  parse_args "$@"

  [[ "$READ_LENGTH" =~ ^[0-9]+$ ]] || die "--read-length must be a positive integer"
  [[ "$THREADS" =~ ^[0-9]+$ ]] || die "--threads must be a positive integer"

  have gzip || die "gzip is required"
  have "$STAR_BIN" || die "STAR not found: $STAR_BIN"

  local species_cap="${SPECIES^}"
  local fasta_url="https://ftp.ensembl.org/pub/release-${RELEASE}/fasta/${SPECIES}/dna/${species_cap}.${ASSEMBLY}.dna.primary_assembly.fa.gz"
  local gtf_url="https://ftp.ensembl.org/pub/release-${RELEASE}/gtf/${SPECIES}/${species_cap}.${ASSEMBLY}.${RELEASE}.gtf.gz"

  local fasta_dir="$OUT_DIR/fasta"
  local gtf_dir="$OUT_DIR/gtf"
  local star_dir="$OUT_DIR/star_index"
  local fasta_gz="$fasta_dir/genome.fa.gz"
  local gtf_gz="$gtf_dir/annotation.gtf.gz"
  local gtf_plain="$OUT_DIR/annotation.gtf"

  mkdir -p "$fasta_dir" "$gtf_dir" "$star_dir"

  echo "[INFO] Downloading FASTA: $fasta_url"
  download_file "$fasta_url" "$fasta_gz"

  echo "[INFO] Downloading GTF: $gtf_url"
  download_file "$gtf_url" "$gtf_gz"

  echo "[INFO] Preparing uncompressed GTF: $gtf_plain"
  gzip -dc "$gtf_gz" > "$gtf_plain"

  local sjdb_overhang=$((READ_LENGTH - 1))
  if [[ "$sjdb_overhang" -lt 1 ]]; then
    die "Invalid --read-length: must be >= 2"
  fi

  echo "[INFO] Building STAR index in: $star_dir"
  "$STAR_BIN" \
    --runThreadN "$THREADS" \
    --runMode genomeGenerate \
    --genomeDir "$star_dir" \
    --genomeFastaFiles "$fasta_gz" \
    --sjdbGTFfile "$gtf_plain" \
    --sjdbOverhang "$sjdb_overhang"

  echo "[OK] Reference prepared"
  echo "[OK] STAR index: $star_dir"
  echo "[OK] GTF: $gtf_plain"
  echo
  echo "Use with pipeline:"
  echo "  --genome-dir $star_dir --gtf $gtf_plain"
}

main "$@"
