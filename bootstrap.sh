#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

ENV_FILE="$ROOT_DIR/.env"
ENV_EXAMPLE="$ROOT_DIR/.env.example"
PRECHECK_SCRIPT="$ROOT_DIR/scripts/preflight.sh"
CONDA_ENV_NAME="star-rseqc"
MODE=""
NON_INTERACTIVE=0
AUTO_RUN=0
RUN_ANALYSIS=0
FASTQ_DIR=""
OUTPUT_DIR=""
GENOME_DIR=""
GTF_FILE=""
JOBS=""
THREADS=""
REF_SETUP_MODE=""
REF_DIR=""
REF_SOURCE_DIR=""
REF_FASTA_FILE=""
REF_READ_LENGTH="101"
ENSEMBL_RELEASE="113"
ENSEMBL_SPECIES="homo_sapiens"
ENSEMBL_ASSEMBLY="GRCh38"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log() { printf "%b[%s]%b %s\n" "$BLUE" "INFO" "$NC" "$*"; }
warn() { printf "%b[%s]%b %s\n" "$YELLOW" "WARN" "$NC" "$*"; }
err() { printf "%b[%s]%b %s\n" "$RED" "ERROR" "$NC" "$*" >&2; }
success() { printf "%b[%s]%b %s\n" "$GREEN" "OK" "$NC" "$*"; }

die() {
  err "$*"
  exit 1
}

have() {
  command -v "$1" >/dev/null 2>&1
}

usage() {
  cat <<USAGE
Usage: ./bootstrap.sh [options]

Options:
  --mode <docker|conda|manual>  Preselect setup mode (skip prompt)
  --non-interactive             Do not prompt; requires --mode
  --run                         Run analysis after setup (dry-run + optional full run)
  --fastq-dir <DIR>             FASTQ input directory (host path)
  --output-dir <DIR>            Output directory (host path)
  --genome-dir <DIR>            STAR genome index directory (host path)
  --gtf <FILE>                  GTF annotation file (host path)
  --jobs <N>                    Parallel jobs override for full run
  --threads <N>                 Threads override for full run
  --ref-mode <existing|local|ensembl>
                                Reference setup mode
  --ref-dir <DIR>               Reference base directory (auto-detects files)
  --ref-source-dir <DIR>        For --ref-mode local: dir containing FASTA + GTF
  --ref-fasta <FILE>            For --ref-mode local: genome FASTA/FA/FNA(.gz)
  --read-length <N>             Read length (for STAR sjdbOverhang, default: 101)
  --ensembl-release <N>         Ensembl release for download mode (default: 113)
  -h, --help                    Show this help
USAGE
}

parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --mode)
        MODE="${2:-}"
        shift 2
        ;;
      --non-interactive)
        NON_INTERACTIVE=1
        shift
        ;;
      --run)
        AUTO_RUN=1
        RUN_ANALYSIS=1
        shift
        ;;
      --fastq-dir)
        FASTQ_DIR="${2:-}"
        RUN_ANALYSIS=1
        shift 2
        ;;
      --output-dir)
        OUTPUT_DIR="${2:-}"
        RUN_ANALYSIS=1
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
      --jobs)
        JOBS="${2:-}"
        shift 2
        ;;
      --threads)
        THREADS="${2:-}"
        shift 2
        ;;
      --ref-mode)
        REF_SETUP_MODE="${2:-}"
        shift 2
        ;;
      --ref-dir)
        REF_DIR="${2:-}"
        shift 2
        ;;
      --ref-source-dir)
        REF_SOURCE_DIR="${2:-}"
        shift 2
        ;;
      --ref-fasta)
        REF_FASTA_FILE="${2:-}"
        shift 2
        ;;
      --read-length)
        REF_READ_LENGTH="${2:-}"
        shift 2
        ;;
      --ensembl-release)
        ENSEMBL_RELEASE="${2:-}"
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

  if [[ $NON_INTERACTIVE -eq 1 && -z "$MODE" ]]; then
    die "--non-interactive requires --mode"
  fi

  if [[ -n "$MODE" ]]; then
    case "$MODE" in
      docker|conda|manual) ;;
      *) die "Invalid --mode '$MODE'. Use docker, conda, or manual." ;;
    esac
  fi

  if [[ -n "$REF_SETUP_MODE" ]]; then
    case "$REF_SETUP_MODE" in
      existing|local|ensembl) ;;
      *) die "Invalid --ref-mode '$REF_SETUP_MODE'. Use existing, local, or ensembl." ;;
    esac
  fi
}

resolve_path() {
  local p="$1"
  if [[ -z "$p" ]]; then
    return 0
  fi
  if [[ "$p" = /* ]]; then
    printf "%s" "$p"
  else
    printf "%s/%s" "$ROOT_DIR" "$p"
  fi
}

detect_state() {
  local docker_ok="no"
  local compose_ok="no"
  local conda_ok="no"
  local mamba_ok="no"
  local cargo_ok="no"

  have docker && docker_ok="yes"
  if have docker && docker compose version >/dev/null 2>&1; then
    compose_ok="yes"
  fi
  have conda && conda_ok="yes"
  have mamba && mamba_ok="yes"
  have cargo && cargo_ok="yes"

  echo
  log "Environment detection"
  printf "  docker:         %s\n" "$docker_ok"
  printf "  docker compose: %s\n" "$compose_ok"
  printf "  conda:          %s\n" "$conda_ok"
  printf "  mamba:          %s\n" "$mamba_ok"
  printf "  cargo:          %s\n" "$cargo_ok"
  echo
}

choose_mode() {
  if [[ -n "$MODE" ]]; then
    return 0
  fi

  if [[ $NON_INTERACTIVE -eq 1 ]]; then
    die "Mode was not provided in non-interactive mode"
  fi

  echo "How do you want to install and run the pipeline?"
  echo "  1) Docker (recommended)"
  echo "  2) Conda"
  echo "  3) Manual install"
  read -r -p "Enter choice [1-3]: " choice
  case "$choice" in
    1) MODE="docker" ;;
    2) MODE="conda" ;;
    3) MODE="manual" ;;
    *) die "Invalid choice: $choice" ;;
  esac
}

prompt_default() {
  local label="$1"
  local default="$2"
  local value
  read -r -p "$label [$default]: " value
  if [[ -z "$value" ]]; then
    printf "%s" "$default"
  else
    printf "%s" "$value"
  fi
}

first_match() {
  local dir="$1"
  shift
  local pattern
  for pattern in "$@"; do
    local f
    f="$(find "$dir" -maxdepth 4 -type f -iname "$pattern" | sort | head -n 1 || true)"
    if [[ -n "$f" ]]; then
      printf "%s" "$f"
      return 0
    fi
  done
  return 1
}

normalize_gtf_if_gz() {
  if [[ -n "$GTF_FILE" && "$GTF_FILE" == *.gz ]]; then
    have gzip || die "gzip is required to unpack GTF: $GTF_FILE"
    local out_gtf="${GTF_FILE%.gz}"
    log "Detected compressed GTF, unpacking: $GTF_FILE -> $out_gtf"
    gzip -dc "$GTF_FILE" > "$out_gtf"
    GTF_FILE="$out_gtf"
    success "Using uncompressed GTF: $GTF_FILE"
  fi
}

interactive_reference_prompt() {
  if [[ -n "$REF_SETUP_MODE" ]]; then
    return 0
  fi

  echo
  echo "Reference genome setup:"
  echo "  1) I already have STAR index + GTF"
  echo "  2) I have local FASTA/FNA/FA + GTF (build STAR index now)"
  echo "  3) Download human reference from Ensembl and build STAR index"
  read -r -p "Enter choice [1-3]: " ref_choice
  case "$ref_choice" in
    1) REF_SETUP_MODE="existing" ;;
    2) REF_SETUP_MODE="local" ;;
    3) REF_SETUP_MODE="ensembl" ;;
    *) die "Invalid reference choice: $ref_choice" ;;
  esac
}

interactive_path_prompt() {
  if [[ $NON_INTERACTIVE -eq 1 ]]; then
    return 0
  fi

  echo
  interactive_reference_prompt

  case "$REF_SETUP_MODE" in
    existing)
      REF_DIR="$(prompt_default "Reference directory (contains STAR index + GTF)" "${REF_DIR:-$ROOT_DIR/refs}")"
      if [[ -d "$REF_DIR" ]]; then
        local auto_gtf auto_genome
        auto_gtf="$(first_match "$REF_DIR" "*.gtf" "*.gtf.gz" || true)"
        auto_genome=""
        if [[ -e "$REF_DIR/Genome" ]]; then
          auto_genome="$REF_DIR"
        elif [[ -e "$REF_DIR/star_index/Genome" ]]; then
          auto_genome="$REF_DIR/star_index"
        else
          local genome_file
          genome_file="$(find "$REF_DIR" -maxdepth 3 -type f -name "Genome" | sort | head -n 1 || true)"
          if [[ -n "$genome_file" ]]; then
            auto_genome="$(dirname "$genome_file")"
          fi
        fi
        if [[ -n "$auto_gtf" ]]; then
          GTF_FILE="$auto_gtf"
        fi
        if [[ -n "$auto_genome" ]]; then
          GENOME_DIR="$auto_genome"
        fi
      fi
      GENOME_DIR="$(prompt_default "STAR genome index directory" "${GENOME_DIR:-$REF_DIR/star_index}")"
      GTF_FILE="$(prompt_default "GTF file path" "${GTF_FILE:-$REF_DIR/annotation.gtf}")"
      ;;
    local)
      REF_DIR="$(prompt_default "Reference directory (contains FASTA/FNA/FA + GTF)" "${REF_DIR:-$ROOT_DIR/refs}")"
      REF_SOURCE_DIR="$REF_DIR"
      if [[ -d "$REF_SOURCE_DIR" ]]; then
        local auto_gtf auto_fa
        auto_gtf="$(first_match "$REF_SOURCE_DIR" "*.gtf" "*.gtf.gz" || true)"
        auto_fa="$(first_match "$REF_SOURCE_DIR" "*.fa" "*.fasta" "*.fna" "*.fa.gz" "*.fasta.gz" "*.fna.gz" || true)"
        if [[ -n "$auto_gtf" ]]; then
          GTF_FILE="$auto_gtf"
        fi
        if [[ -n "$auto_fa" ]]; then
          REF_FASTA_FILE="$auto_fa"
        fi
      fi
      GTF_FILE="$(prompt_default "GTF file path" "${GTF_FILE:-$REF_DIR/annotation.gtf}")"
      REF_FASTA_FILE="$(prompt_default "FASTA/FNA/FA file path" "${REF_FASTA_FILE:-$REF_DIR/genome.fa.gz}")"
      GENOME_DIR="$(prompt_default "STAR index output directory" "${GENOME_DIR:-$REF_DIR/star_index}")"
      REF_READ_LENGTH="$(prompt_default "Read length (for sjdbOverhang)" "$REF_READ_LENGTH")"
      ;;
    ensembl)
      REF_DIR="$(prompt_default "Reference output base directory" "${REF_DIR:-$ROOT_DIR/refs}")"
      GENOME_DIR="$REF_DIR/star_index"
      GTF_FILE="$REF_DIR/annotation.gtf"
      ENSEMBL_RELEASE="$(prompt_default "Ensembl release" "$ENSEMBL_RELEASE")"
      REF_READ_LENGTH="$(prompt_default "Read length (for sjdbOverhang)" "$REF_READ_LENGTH")"
      ;;
    *)
      die "Unexpected REF_SETUP_MODE: $REF_SETUP_MODE"
      ;;
  esac
}

ensure_analysis_paths() {
  if [[ $RUN_ANALYSIS -ne 1 ]]; then
    return 0
  fi

  if [[ -z "$OUTPUT_DIR" ]]; then
    OUTPUT_DIR="$ROOT_DIR/results"
  fi

  if [[ -z "$FASTQ_DIR" ]]; then
    if [[ $NON_INTERACTIVE -eq 1 ]]; then
      die "--fastq-dir is required when --run/analysis is requested in non-interactive mode"
    fi
    FASTQ_DIR="$(prompt_default "FASTQ directory (contains *_1P.fastq.gz, *_2P.fastq.gz)" "$ROOT_DIR/data")"
  fi
  if [[ $NON_INTERACTIVE -eq 0 ]]; then
    OUTPUT_DIR="$(prompt_default "Output directory" "$OUTPUT_DIR")"
  fi

  FASTQ_DIR="$(resolve_path "$FASTQ_DIR")"
  OUTPUT_DIR="$(resolve_path "$OUTPUT_DIR")"
}

build_star_index_local() {
  local star_bin="$1"
  [[ -f "$GTF_FILE" ]] || die "GTF file not found: $GTF_FILE"
  [[ -f "$REF_FASTA_FILE" ]] || die "FASTA file not found: $REF_FASTA_FILE"
  [[ "$REF_READ_LENGTH" =~ ^[0-9]+$ ]] || die "--read-length must be a positive integer"
  local sjdb=$((REF_READ_LENGTH - 1))
  [[ "$sjdb" -ge 1 ]] || die "Read length must be >= 2"

  mkdir -p "$GENOME_DIR"
  log "Building STAR index from local FASTA/GTF"
  "$star_bin" \
    --runThreadN "$(nproc 2>/dev/null || echo 8)" \
    --runMode genomeGenerate \
    --genomeDir "$GENOME_DIR" \
    --genomeFastaFiles "$REF_FASTA_FILE" \
    --sjdbGTFfile "$GTF_FILE" \
    --sjdbOverhang "$sjdb"
  success "STAR index ready: $GENOME_DIR"
}

prepare_references_host() {
  case "${REF_SETUP_MODE:-existing}" in
    existing)
      return 0
      ;;
    local)
      have STAR || die "STAR is required to build a local index. Install STAR or use --ref-mode existing."
      build_star_index_local "STAR"
      ;;
    ensembl)
      have STAR || die "STAR is required to download/index Ensembl references. Install STAR or use --ref-mode existing."
      local out_base
      out_base="$(dirname "$GTF_FILE")"
      log "Downloading Ensembl reference (release ${ENSEMBL_RELEASE}) and building STAR index"
      "$ROOT_DIR/scripts/download_and_index_reference.sh" \
        --release "$ENSEMBL_RELEASE" \
        --assembly "$ENSEMBL_ASSEMBLY" \
        --species "$ENSEMBL_SPECIES" \
        --read-length "$REF_READ_LENGTH" \
        --out-dir "$out_base" \
        --star-bin STAR
      GENOME_DIR="$out_base/star_index"
      GTF_FILE="$out_base/annotation.gtf"
      success "Reference prepared in: $out_base"
      ;;
    *)
      die "Unexpected REF_SETUP_MODE: $REF_SETUP_MODE"
      ;;
  esac
}

prepare_references_conda() {
  local ccmd="$1"
  case "${REF_SETUP_MODE:-existing}" in
    existing)
      return 0
      ;;
    local)
      [[ -f "$GTF_FILE" ]] || die "GTF file not found: $GTF_FILE"
      [[ -f "$REF_FASTA_FILE" ]] || die "FASTA file not found: $REF_FASTA_FILE"
      [[ "$REF_READ_LENGTH" =~ ^[0-9]+$ ]] || die "--read-length must be a positive integer"
      local sjdb=$((REF_READ_LENGTH - 1))
      [[ "$sjdb" -ge 1 ]] || die "Read length must be >= 2"
      mkdir -p "$GENOME_DIR"
      log "Building STAR index in conda environment"
      "$ccmd" run -n "$CONDA_ENV_NAME" STAR \
        --runThreadN "$(nproc 2>/dev/null || echo 8)" \
        --runMode genomeGenerate \
        --genomeDir "$GENOME_DIR" \
        --genomeFastaFiles "$REF_FASTA_FILE" \
        --sjdbGTFfile "$GTF_FILE" \
        --sjdbOverhang "$sjdb"
      success "STAR index ready: $GENOME_DIR"
      ;;
    ensembl)
      local out_base
      out_base="$(dirname "$GTF_FILE")"
      log "Downloading Ensembl reference inside conda environment"
      "$ccmd" run -n "$CONDA_ENV_NAME" "$ROOT_DIR/scripts/download_and_index_reference.sh" \
        --release "$ENSEMBL_RELEASE" \
        --assembly "$ENSEMBL_ASSEMBLY" \
        --species "$ENSEMBL_SPECIES" \
        --read-length "$REF_READ_LENGTH" \
        --out-dir "$out_base" \
        --star-bin STAR
      GENOME_DIR="$out_base/star_index"
      GTF_FILE="$out_base/annotation.gtf"
      success "Reference prepared in: $out_base"
      ;;
    *)
      die "Unexpected REF_SETUP_MODE: $REF_SETUP_MODE"
      ;;
  esac
}

prepare_references_docker() {
  case "${REF_SETUP_MODE:-existing}" in
    existing)
      return 0
      ;;
    local)
      [[ -f "$GTF_FILE" ]] || die "GTF file not found: $GTF_FILE"
      [[ -f "$REF_FASTA_FILE" ]] || die "FASTA file not found: $REF_FASTA_FILE"
      [[ "$REF_READ_LENGTH" =~ ^[0-9]+$ ]] || die "--read-length must be a positive integer"
      local sjdb=$((REF_READ_LENGTH - 1))
      [[ "$sjdb" -ge 1 ]] || die "Read length must be >= 2"

      local gtf_parent gtf_name fasta_parent fasta_name genome_parent genome_name threads
      gtf_parent="$(dirname "$GTF_FILE")"
      gtf_name="$(basename "$GTF_FILE")"
      fasta_parent="$(dirname "$REF_FASTA_FILE")"
      fasta_name="$(basename "$REF_FASTA_FILE")"
      genome_parent="$(dirname "$GENOME_DIR")"
      genome_name="$(basename "$GENOME_DIR")"
      threads="$(nproc 2>/dev/null || echo 8)"
      mkdir -p "$genome_parent"

      log "Building STAR index in Docker container"
      docker compose run --rm \
        --entrypoint STAR \
        -v "$gtf_parent:/gtfsrc:ro" \
        -v "$fasta_parent:/fastasrc:ro" \
        -v "$genome_parent:/outroot" \
        star-rseqc \
        --runThreadN "$threads" \
        --runMode genomeGenerate \
        --genomeDir "/outroot/$genome_name" \
        --genomeFastaFiles "/fastasrc/$fasta_name" \
        --sjdbGTFfile "/gtfsrc/$gtf_name" \
        --sjdbOverhang "$sjdb"
      success "STAR index ready: $GENOME_DIR"
      ;;
    ensembl)
      local out_base
      out_base="$(dirname "$GTF_FILE")"
      mkdir -p "$out_base"
      log "Downloading Ensembl reference and building STAR index in Docker container"
      docker compose run --rm \
        --entrypoint /bin/sh \
        -v "$ROOT_DIR/scripts:/scripts:ro" \
        -v "$out_base:/refsout" \
        star-rseqc \
        -lc "/scripts/download_and_index_reference.sh \
          --release '$ENSEMBL_RELEASE' \
          --assembly '$ENSEMBL_ASSEMBLY' \
          --species '$ENSEMBL_SPECIES' \
          --read-length '$REF_READ_LENGTH' \
          --out-dir /refsout \
          --star-bin STAR"
      GENOME_DIR="$out_base/star_index"
      GTF_FILE="$out_base/annotation.gtf"
      success "Reference prepared in: $out_base"
      ;;
    *)
      die "Unexpected REF_SETUP_MODE: $REF_SETUP_MODE"
      ;;
  esac
}

apply_ref_dir_hints() {
  if [[ -z "$REF_DIR" || ! -d "$REF_DIR" ]]; then
    return 0
  fi

  case "${REF_SETUP_MODE:-existing}" in
    existing)
      if [[ -z "$GENOME_DIR" ]]; then
        if [[ -e "$REF_DIR/Genome" ]]; then
          GENOME_DIR="$REF_DIR"
        elif [[ -e "$REF_DIR/star_index/Genome" ]]; then
          GENOME_DIR="$REF_DIR/star_index"
        else
          local genome_file
          genome_file="$(find "$REF_DIR" -maxdepth 3 -type f -name "Genome" | sort | head -n 1 || true)"
          if [[ -n "$genome_file" ]]; then
            GENOME_DIR="$(dirname "$genome_file")"
          fi
        fi
      fi
      if [[ -z "$GTF_FILE" ]]; then
        GTF_FILE="$(first_match "$REF_DIR" "*.gtf" "*.gtf.gz" || true)"
      fi
      ;;
    local)
      if [[ -z "$REF_SOURCE_DIR" ]]; then
        REF_SOURCE_DIR="$REF_DIR"
      fi
      if [[ -z "$GTF_FILE" ]]; then
        GTF_FILE="$(first_match "$REF_DIR" "*.gtf" "*.gtf.gz" || true)"
      fi
      if [[ -z "$REF_FASTA_FILE" ]]; then
        REF_FASTA_FILE="$(first_match "$REF_DIR" "*.fa" "*.fasta" "*.fna" "*.fa.gz" "*.fasta.gz" "*.fna.gz" || true)"
      fi
      if [[ -z "$GENOME_DIR" ]]; then
        GENOME_DIR="$REF_DIR/star_index"
      fi
      ;;
    ensembl)
      if [[ -z "$GENOME_DIR" ]]; then
        GENOME_DIR="$REF_DIR/star_index"
      fi
      if [[ -z "$GTF_FILE" ]]; then
        GTF_FILE="$REF_DIR/annotation.gtf"
      fi
      ;;
    *)
      ;;
  esac
}

should_run_full() {
  if [[ $AUTO_RUN -eq 1 ]]; then
    return 0
  fi

  if [[ $NON_INTERACTIVE -eq 1 ]]; then
    return 1
  fi

  read -r -p "Dry-run passed. Start full analysis now? [y/N]: " answer
  case "$answer" in
    y|Y|yes|YES) return 0 ;;
    *) return 1 ;;
  esac
}

ensure_precheck_script() {
  [[ -x "$PRECHECK_SCRIPT" ]] || die "Missing executable preflight script: $PRECHECK_SCRIPT"
}

load_path_defaults() {
  # Load optional .env defaults while preserving CLI-provided values.
  if [[ -f "$ENV_FILE" ]]; then
    local cli_fastq="$FASTQ_DIR"
    local cli_output="$OUTPUT_DIR"
    local cli_genome="$GENOME_DIR"
    local cli_gtf="$GTF_FILE"
    local cli_ref_dir="$REF_DIR"
    local cli_jobs="$JOBS"
    local cli_threads="$THREADS"
    while IFS= read -r line || [[ -n "$line" ]]; do
      # trim leading/trailing spaces
      line="${line#"${line%%[![:space:]]*}"}"
      line="${line%"${line##*[![:space:]]}"}"
      [[ -z "$line" || "$line" == \#* ]] && continue
      [[ "$line" == *=* ]] || continue

      local key="${line%%=*}"
      local val="${line#*=}"
      key="${key#"${key%%[![:space:]]*}"}"
      key="${key%"${key##*[![:space:]]}"}"
      val="${val#"${val%%[![:space:]]*}"}"
      val="${val%"${val##*[![:space:]]}"}"

      # Strip single/double quotes if present around the full value.
      if [[ "$val" == \"*\" && "$val" == *\" ]]; then
        val="${val:1:${#val}-2}"
      elif [[ "$val" == \'*\' && "$val" == *\' ]]; then
        val="${val:1:${#val}-2}"
      fi

      case "$key" in
        FASTQ_DIR|OUTPUT_DIR|GENOME_DIR|GTF_FILE|REF_DIR|JOBS|THREADS)
          printf -v "$key" '%s' "$val"
          ;;
        *)
          ;;
      esac
    done < "$ENV_FILE"

    if [[ -n "$cli_fastq" ]]; then FASTQ_DIR="$cli_fastq"; fi
    if [[ -n "$cli_output" ]]; then OUTPUT_DIR="$cli_output"; fi
    if [[ -n "$cli_genome" ]]; then GENOME_DIR="$cli_genome"; fi
    if [[ -n "$cli_gtf" ]]; then GTF_FILE="$cli_gtf"; fi
    if [[ -n "$cli_ref_dir" ]]; then REF_DIR="$cli_ref_dir"; fi
    if [[ -n "$cli_jobs" ]]; then JOBS="$cli_jobs"; fi
    if [[ -n "$cli_threads" ]]; then THREADS="$cli_threads"; fi
  fi

  if [[ -z "$GENOME_DIR" ]]; then
    GENOME_DIR="$ROOT_DIR/refs/star_index"
  fi
  if [[ -z "$GTF_FILE" ]]; then
    GTF_FILE="$ROOT_DIR/refs/annotation.gtf"
  fi

  # Allow .env values using container-style /work/* paths.
  if [[ "$FASTQ_DIR" == /work/* ]]; then FASTQ_DIR="$ROOT_DIR/${FASTQ_DIR#/work/}"; fi
  if [[ "$OUTPUT_DIR" == /work/* ]]; then OUTPUT_DIR="$ROOT_DIR/${OUTPUT_DIR#/work/}"; fi
  if [[ "$GENOME_DIR" == /work/* ]]; then GENOME_DIR="$ROOT_DIR/${GENOME_DIR#/work/}"; fi
  if [[ "$GTF_FILE" == /work/* ]]; then GTF_FILE="$ROOT_DIR/${GTF_FILE#/work/}"; fi

  if [[ -n "$FASTQ_DIR" ]]; then
    FASTQ_DIR="$(resolve_path "$FASTQ_DIR")"
  fi
  if [[ -n "$OUTPUT_DIR" ]]; then
    OUTPUT_DIR="$(resolve_path "$OUTPUT_DIR")"
  fi
  GENOME_DIR="$(resolve_path "$GENOME_DIR")"
  GTF_FILE="$(resolve_path "$GTF_FILE")"
  if [[ -n "$REF_SOURCE_DIR" ]]; then
    REF_SOURCE_DIR="$(resolve_path "$REF_SOURCE_DIR")"
  fi
  if [[ -n "$REF_DIR" ]]; then
    REF_DIR="$(resolve_path "$REF_DIR")"
  fi
  if [[ -n "$REF_FASTA_FILE" ]]; then
    REF_FASTA_FILE="$(resolve_path "$REF_FASTA_FILE")"
  fi
}

setup_repo_dirs() {
  mkdir -p "$ROOT_DIR/data" "$ROOT_DIR/refs" "$ROOT_DIR/results" "$ROOT_DIR/scripts"
}

ensure_env_file() {
  if [[ ! -f "$ENV_FILE" && -f "$ENV_EXAMPLE" ]]; then
    cp "$ENV_EXAMPLE" "$ENV_FILE"
    success "Created .env from .env.example"
  fi
}

run_preflight() {
  "$PRECHECK_SCRIPT" \
    --fastq-dir "$FASTQ_DIR" \
    --output-dir "$OUTPUT_DIR" \
    --genome-dir "$GENOME_DIR" \
    --gtf "$GTF_FILE"
}

append_runtime_flags() {
  local -n _arr_ref=$1
  if [[ -n "$JOBS" ]]; then
    _arr_ref+=(--jobs "$JOBS")
  fi
  if [[ -n "$THREADS" ]]; then
    _arr_ref+=(--threads "$THREADS")
  fi
}

run_docker_mode() {
  have docker || die "Docker is not installed. Install Docker Desktop/Engine and re-run bootstrap."
  docker compose version >/dev/null 2>&1 || die "Docker Compose plugin is missing. Install it and re-run bootstrap."

  setup_repo_dirs
  ensure_env_file
  load_path_defaults
  apply_ref_dir_hints
  interactive_path_prompt
  log "Building container image"
  docker compose build

  prepare_references_docker
  normalize_gtf_if_gz
  if [[ $RUN_ANALYSIS -ne 1 ]]; then
    success "Setup complete. Dependencies and references are ready."
    log "Run analysis later with: ./setup.sh --mode docker --run --fastq-dir ./data --output-dir ./results"
    return 0
  fi
  ensure_analysis_paths
  ensure_precheck_script

  log "Running host preflight checks"
  run_preflight

  local gtf_parent
  local gtf_name
  gtf_parent="$(dirname "$GTF_FILE")"
  gtf_name="$(basename "$GTF_FILE")"

  local docker_args=(--rm
    -v "$FASTQ_DIR:/input:ro"
    -v "$GENOME_DIR:/genome:ro"
    -v "$gtf_parent:/gtf:ro"
    -v "$OUTPUT_DIR:/output")

  local dry_args=(/input --output /output --genome-dir /genome --gtf "/gtf/$gtf_name" --dry-run)
  append_runtime_flags dry_args

  log "Running pipeline dry-run in container"
  docker compose run "${docker_args[@]}" star-rseqc "${dry_args[@]}"

  if should_run_full; then
    local run_args=(/input --output /output --genome-dir /genome --gtf "/gtf/$gtf_name")
    append_runtime_flags run_args
    log "Starting full analysis in container"
    docker compose run "${docker_args[@]}" star-rseqc "${run_args[@]}"
  else
    log "Skipping full run. Re-run with: ./bootstrap.sh --mode docker --run"
  fi
}

choose_conda_cmd() {
  if have mamba; then
    echo "mamba"
  elif have conda; then
    echo "conda"
  else
    echo ""
  fi
}

conda_env_exists() {
  local ccmd="$1"
  "$ccmd" env list | awk '{print $1}' | grep -qx "$CONDA_ENV_NAME"
}

run_conda_mode() {
  local ccmd
  ccmd="$(choose_conda_cmd)"
  [[ -n "$ccmd" ]] || die "Conda or Mamba is required for conda mode."

  setup_repo_dirs
  ensure_env_file
  load_path_defaults
  apply_ref_dir_hints
  interactive_path_prompt
  ensure_precheck_script

  log "Preparing conda environment: $CONDA_ENV_NAME"
  if conda_env_exists "$ccmd"; then
    "$ccmd" env update -n "$CONDA_ENV_NAME" -f "$ROOT_DIR/environment.yml" --prune
  else
    "$ccmd" env create -n "$CONDA_ENV_NAME" -f "$ROOT_DIR/environment.yml"
  fi

  prepare_references_conda "$ccmd"
  normalize_gtf_if_gz
  if [[ $RUN_ANALYSIS -ne 1 ]]; then
    success "Setup complete. Conda environment and references are ready."
    log "Run analysis later with: ./setup.sh --mode conda --run --fastq-dir ./data --output-dir ./results"
    return 0
  fi
  ensure_analysis_paths

  log "Running host preflight checks"
  run_preflight

  log "Building Rust binary in conda environment"
  "$ccmd" run -n "$CONDA_ENV_NAME" cargo build --release

  local run_cmd=("$ROOT_DIR/target/release/star-rseqc" "$FASTQ_DIR" --output "$OUTPUT_DIR" --genome-dir "$GENOME_DIR" --gtf "$GTF_FILE" --dry-run)
  append_runtime_flags run_cmd

  log "Running pipeline dry-run in conda environment"
  local cprefix
  cprefix="$("$ccmd" run -n "$CONDA_ENV_NAME" python -c 'import os; print(os.environ["CONDA_PREFIX"])')"
  "$ccmd" run -n "$CONDA_ENV_NAME" env \
    STAR_RSEQC_STAR_ENV="$cprefix" \
    STAR_RSEQC_RSEQC_ENV="$cprefix" \
    STAR_RSEQC_DEEPTOOLS_ENV="$cprefix" \
    "${run_cmd[@]}"

  if should_run_full; then
    local full_cmd=("$ROOT_DIR/target/release/star-rseqc" "$FASTQ_DIR" --output "$OUTPUT_DIR" --genome-dir "$GENOME_DIR" --gtf "$GTF_FILE")
    append_runtime_flags full_cmd
    log "Starting full analysis in conda environment"
    "$ccmd" run -n "$CONDA_ENV_NAME" env \
      STAR_RSEQC_STAR_ENV="$cprefix" \
      STAR_RSEQC_RSEQC_ENV="$cprefix" \
      STAR_RSEQC_DEEPTOOLS_ENV="$cprefix" \
      "${full_cmd[@]}"
  else
    log "Skipping full run. Re-run with --mode conda --run when ready."
  fi
}

run_manual_mode() {
  setup_repo_dirs
  ensure_env_file
  load_path_defaults
  apply_ref_dir_hints
  interactive_path_prompt
  prepare_references_host
  normalize_gtf_if_gz
  if [[ $RUN_ANALYSIS -ne 1 ]]; then
    success "Setup complete. References are ready."
    log "Run analysis later with: ./setup.sh --mode manual --run --fastq-dir ./data --output-dir ./results"
    return 0
  fi
  ensure_analysis_paths
  ensure_precheck_script

  local missing=0
  local tools=(STAR samtools infer_experiment.py read_distribution.py geneBody_coverage2.py bamCoverage Rscript cargo)

  log "Checking required tools for manual mode"
  for t in "${tools[@]}"; do
    if have "$t"; then
      printf "  %-24s %s\n" "$t" "found"
    else
      printf "  %-24s %s\n" "$t" "missing"
      missing=1
    fi
  done

  if have Rscript; then
    if ! Rscript -e 'library(ggplot2)' >/dev/null 2>&1; then
      warn "R package ggplot2 is missing"
      missing=1
    fi
  fi

  if [[ $missing -ne 0 ]]; then
    cat <<'MANUAL_HELP'

Manual mode prerequisites are incomplete.
Install missing tools, then re-run bootstrap.

Ubuntu/Debian base packages:
  sudo apt-get update
  sudo apt-get install -y build-essential curl git samtools r-base

R package:
  Rscript -e 'install.packages("ggplot2", repos="https://cloud.r-project.org")'

Recommended for full bio stack:
  Use --mode conda or --mode docker for automated STAR/RSeQC/deeptools setup.
MANUAL_HELP
    exit 2
  fi

  log "Running host preflight checks"
  run_preflight

  log "Building Rust binary"
  cargo build --release

  local dry_cmd=("$ROOT_DIR/target/release/star-rseqc" "$FASTQ_DIR" --output "$OUTPUT_DIR" --genome-dir "$GENOME_DIR" --gtf "$GTF_FILE" --dry-run)
  append_runtime_flags dry_cmd

  log "Running pipeline dry-run"
  "${dry_cmd[@]}"

  if should_run_full; then
    local full_cmd=("$ROOT_DIR/target/release/star-rseqc" "$FASTQ_DIR" --output "$OUTPUT_DIR" --genome-dir "$GENOME_DIR" --gtf "$GTF_FILE")
    append_runtime_flags full_cmd
    log "Starting full analysis"
    "${full_cmd[@]}"
  else
    log "Skipping full run. Launch later with target/release/star-rseqc ..."
  fi
}

main() {
  parse_args "$@"
  detect_state
  choose_mode

  case "$MODE" in
    docker) run_docker_mode ;;
    conda) run_conda_mode ;;
    manual) run_manual_mode ;;
    *) die "Unexpected mode: $MODE" ;;
  esac

  success "Bootstrap completed"
}

main "$@"
