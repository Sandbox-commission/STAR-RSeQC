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
FASTQ_DIR=""
OUTPUT_DIR=""
GENOME_DIR=""
GTF_FILE=""
JOBS=""
THREADS=""

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
  --run                         Automatically run full analysis after dry-run
  --fastq-dir <DIR>             FASTQ input directory (host path)
  --output-dir <DIR>            Output directory (host path)
  --genome-dir <DIR>            STAR genome index directory (host path)
  --gtf <FILE>                  GTF annotation file (host path)
  --jobs <N>                    Parallel jobs override for full run
  --threads <N>                 Threads override for full run
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
        shift
        ;;
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
      --jobs)
        JOBS="${2:-}"
        shift 2
        ;;
      --threads)
        THREADS="${2:-}"
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

  echo "Choose setup mode:"
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
        FASTQ_DIR|OUTPUT_DIR|GENOME_DIR|GTF_FILE|JOBS|THREADS)
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
    if [[ -n "$cli_jobs" ]]; then JOBS="$cli_jobs"; fi
    if [[ -n "$cli_threads" ]]; then THREADS="$cli_threads"; fi
  fi

  if [[ -z "$FASTQ_DIR" ]]; then
    FASTQ_DIR="$ROOT_DIR/data"
  fi
  if [[ -z "$OUTPUT_DIR" ]]; then
    OUTPUT_DIR="$ROOT_DIR/results"
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

  FASTQ_DIR="$(resolve_path "$FASTQ_DIR")"
  OUTPUT_DIR="$(resolve_path "$OUTPUT_DIR")"
  GENOME_DIR="$(resolve_path "$GENOME_DIR")"
  GTF_FILE="$(resolve_path "$GTF_FILE")"
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
  ensure_precheck_script

  log "Running host preflight checks"
  run_preflight

  log "Building container image"
  docker compose build

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
  ensure_precheck_script

  log "Running host preflight checks"
  run_preflight

  log "Preparing conda environment: $CONDA_ENV_NAME"
  if conda_env_exists "$ccmd"; then
    "$ccmd" env update -n "$CONDA_ENV_NAME" -f "$ROOT_DIR/environment.yml" --prune
  else
    "$ccmd" env create -n "$CONDA_ENV_NAME" -f "$ROOT_DIR/environment.yml"
  fi

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
