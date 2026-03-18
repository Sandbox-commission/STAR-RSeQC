#!/usr/bin/env bash
set -euo pipefail

if [[ $# -eq 0 ]]; then
  echo "Usage: ./scripts/run-container.sh <FASTQ_DIR> [star-rseqc options...]" >&2
  echo "Example: ./scripts/run-container.sh /work/data --output /work/results --genome-dir /work/refs/star_index --gtf /work/refs/annotation.gtf" >&2
  exit 1
fi

docker compose run --rm star-rseqc "$@"
