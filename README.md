# star-rseqc

**STAR 2-pass alignment + deeptools + RSeQC quality control pipeline for paired-end RNA-seq**

A high-performance, resume-aware pipeline written in Rust that automates STAR
two-pass alignment, BAM-to-bigwig conversion, and RSeQC quality control with
ggplot2 PDF generation for bulk paired-end RNA-seq experiments. Features a
full-screen terminal UI with real-time progress tracking and per-step SHA256
checkpoint verification.

---

## Table of Contents

- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
  - [Arguments](#arguments)
  - [Options](#options)
- [FASTQ Naming Convention](#fastq-naming-convention)
- [Pipeline Phases](#pipeline-phases)
- [Output Structure](#output-structure)
- [STAR Parameters](#star-parameters)
- [Resume and SHA256 Integrity](#resume-and-sha256-integrity)
  - [How It Works](#how-it-works)
  - [Checkpoint Format](#checkpoint-format)
  - [Resume States](#resume-states)
- [Terminal UI](#terminal-ui)
- [Reference Configuration](#reference-configuration)
- [Examples](#examples)
- [Architecture](#architecture)
- [License](#license)

---

## Features

- **Three-phase pipeline**: STAR alignment → deeptools bigwig → RSeQC QC
- **STAR 2-pass alignment** with chimeric junction detection, transcriptome BAM,
  and gene-level quantification (ENCODE-compliant parameters)
- **deeptools bamCoverage** — BAM-to-bigwig conversion (binSize 10)
- **RSeQC quality control**: strandedness inference, gene body coverage, and
  read distribution analysis — all three tools run in parallel per sample
- **ggplot2 PDF plots** — gene body coverage curves and heatmap PDFs generated
  via Rscript (mandatory, included in checkpoint verification)
- **Pure-Rust GTF-to-BED12 conversion** — no external tools needed; the
  annotation is converted automatically on first run and cached
- **Full-screen TUI** — real-time progress monitor with per-sample spinners,
  overall progress bar, active job slots, elapsed/ETA timers, and a scrolling
  activity log (built with `crossterm`)
- **Per-step SHA256 checkpoints** — each pipeline step (STAR, deeptools, infer,
  read_distribution, geneBody_coverage) writes its own SHA256 checkpoint
  immediately after completion; on resume, only steps with missing or mismatched
  checkpoints are re-run, and partial outputs are cleaned automatically
- **Parallel execution** — configurable per-phase job counts with auto-detection
  from system RAM and CPU count
- **Graceful cancellation** — `Ctrl+C` signals all running jobs to stop cleanly;
  partial outputs are removed so corrupted files never persist
- **Crash-safe checkpoints** — atomic writes with `fsync` ensure checkpoint
  files are never half-written on power loss
- **Dry-run mode** — list discovered samples and their resume status without
  executing anything

---

## Requirements

| Tool | Version | Default Path |
|------|---------|--------------|
| [STAR](https://github.com/alexdobin/STAR) | v2.7.11b+ | `/home/cml/miniforge3/envs/star/bin/STAR` |
| [samtools](http://www.htslib.org/) | v1.15+ | System PATH |
| [RSeQC](http://rseqc.sourceforge.net/) | v5.0+ | `/home/cml/miniforge3/envs/RSeQC/bin/` |
| [deeptools](https://deeptools.readthedocs.io/) | v3.5+ | `/home/cml/miniforge3/envs/deeptools/bin/bamCoverage` |
| [Rscript](https://www.r-project.org/) + ggplot2 | R 4.0+ | System PATH |
| [zcat](https://www.gnu.org/software/gzip/) | any | System PATH |
| Rust toolchain | 1.70+ (edition 2021) | — |

### Reference files

| File | Default Path |
|------|-------------|
| STAR genome index | `/home/cml/humandb/transcriptomeindex/ensembl113/star_hg38_101bp_index` |
| GTF annotation | `/home/cml/humandb/transcriptomeindex/ensembl113/Homo_sapiens.GRCh38.113.gtf` |

All default paths can be overridden via command-line flags or environment variables
(`STAR_RSEQC_GENOME_DIR`, `STAR_RSEQC_GTF`, `STAR_RSEQC_STAR_ENV`,
`STAR_RSEQC_RSEQC_ENV`, `STAR_RSEQC_DEEPTOOLS_ENV`).

---

## Installation

```bash
cd /home/cml/rust-codes/star-rseqc
cargo build --release

# The binary is at:
#   target/release/star-rseqc

# Optionally copy to a directory in your PATH:
cp target/release/star-rseqc ~/.local/bin/
```

### Dependencies (Cargo.toml)

| Crate | Purpose |
|-------|---------|
| `chrono` | Timestamps in logs and summary |
| `crossterm` | Full-screen TUI rendering (alternate screen, raw mode, colors) |
| `glob` | FASTQ file pattern matching |
| `serde` + `serde_json` | JSON summary output |
| `sha2` | SHA256 digest computation for output integrity |

---

## Quick Start

```bash
# Run on a directory containing paired-end FASTQs
star-rseqc /path/to/Paired/

# Run on the current directory
star-rseqc ./

# Dry run to verify sample discovery
star-rseqc ./ --dry-run

# Resume after interruption (just re-run the same command)
star-rseqc ./
```

---

## Usage

```
star-rseqc <FASTQ_DIR> [OPTIONS]
```

### Arguments

| Argument | Description |
|----------|-------------|
| `<FASTQ_DIR>` | Directory containing `*_1P.fastq.gz` paired-end FASTQ files |

### Options

| Flag | Description | Default |
|------|-------------|---------|
| `-o, --output <DIR>` | Output directory | `star-rseqc-results` |
| `-j, --jobs <N>` | Default parallel jobs | auto-detected from RAM |
| `--star-jobs <N>` | STAR phase parallel jobs | same as `--jobs` |
| `--deeptools-jobs <N>` | deeptools phase parallel jobs | same as `--jobs` |
| `--rseqc-jobs <N>` | RSeQC phase parallel jobs | same as `--jobs` |
| `-t, --threads <N>` | Threads per job | auto-detected from CPUs |
| `--bam-sort-ram <BYTES>` | RAM for STAR BAM sorting | auto-detected |
| `--genome-dir <DIR>` | STAR genome index directory | see defaults |
| `--gtf <FILE>` | GTF annotation file | see defaults |
| `--bed <FILE>` | Pre-built BED12 (else auto-generated from GTF) | auto |
| `--star-env <DIR>` | STAR conda environment prefix | see defaults |
| `--rseqc-env <DIR>` | RSeQC conda environment prefix | see defaults |
| `--deeptools-env <DIR>` | deeptools conda environment prefix | see defaults |
| `--samtools <PATH>` | Path to samtools binary | `samtools` |
| `--skip-alignment` | Skip STAR (run QC on existing BAMs) | off |
| `--dry-run` | Preview samples without running | off |
| `--clean-sam` | Delete chimeric SAM files after successful run | off |
| `--operator <NAME>` | Operator name for audit trail | `$USER` |
| `-h, --help` | Print help | — |

---

## FASTQ Naming Convention

Files must follow this pattern:

```
<SAMPLE>_1P.fastq.gz    (read 1 / forward)
<SAMPLE>_2P.fastq.gz    (read 2 / reverse)
```

The sample name is everything before `_1P` or `_2P`:

```
103N_GBC_1P.fastq.gz   →  sample = 103N_GBC
50T_CRC_1P.fastq.gz    →  sample = 50T_CRC
```

Both R1 and R2 must exist for a sample to be included. Sample names are
validated to contain only `[A-Za-z0-9_.-]`.

---

## Pipeline Phases

### Phase 1: STAR 2-Pass Alignment

Runs STAR in two-pass mode with chimeric junction detection:

- Coordinate-sorted BAM (`*_Aligned.sortedByCoord.out.bam`) + BAM index
- Transcriptome BAM (`*_Aligned.toTranscriptome.out.bam`)
- Gene-level counts (`*_ReadsPerGene.out.tab`)
- Chimeric junctions (`*_Chimeric.out.junction`)
- Splice junctions (`*_SJ.out.tab`)
- STAR logs redirected to `logs/<sample>.star.log`

If STAR fails, all partial output files (including `_STARtmp*` directories) are
cleaned up automatically.

### Phase 2: deeptools bamCoverage

Converts each coordinate-sorted BAM to a bigwig file for downstream
visualization and geneBody_coverage2:

- `bigwig/<sample>.bw` (binSize 10, multi-threaded)

### Phase 3: RSeQC Quality Control

Three RSeQC tools run **in parallel** per sample. Each tool writes its own
SHA256 checkpoint on success, so a failure in one tool does not force re-running
the others on resume.

| Tool | Output | Purpose |
|------|--------|---------|
| `infer_experiment.py` | `<sample>.strand.txt` | Library strandedness detection |
| `read_distribution.py` | `<sample>.read_distribution.txt` | Genomic feature distribution |
| `geneBody_coverage2.py` | `<sample>.geneBodyCoverage.{txt,pdf,...}` | 5'-to-3' coverage uniformity |

After geneBody_coverage2 completes, **ggplot2 PDF plots** are generated:

| Plot | File | Description |
|------|------|-------------|
| Curves | `<sample>.geneBodyCoverage.curves.pdf` | Coverage distribution (geom_point + geom_line) |
| Heatmap | `<sample>.geneBodyCoverage.heatMap.pdf` | Normalized coverage intensity (geom_tile) |

PDF generation is mandatory and included in the genebody SHA256 checkpoint.

---

## Output Structure

```
<output>/
├── star/                                    STAR alignment outputs
│   ├── <sample>_Aligned.sortedByCoord.out.bam
│   ├── <sample>_Aligned.sortedByCoord.out.bam.bai
│   ├── <sample>_Aligned.toTranscriptome.out.bam
│   ├── <sample>_ReadsPerGene.out.tab
│   ├── <sample>_SJ.out.tab
│   ├── <sample>_Chimeric.out.junction
│   ├── <sample>_Chimeric.out.sam
│   ├── <sample>_Log.out
│   ├── <sample>_Log.progress.out
│   └── <sample>_Log.final.out
├── bigwig/                                  deeptools bigwig files
│   └── <sample>.bw
├── qc/                                      RSeQC quality control outputs
│   ├── <sample>.strand.txt
│   ├── <sample>.read_distribution.txt
│   ├── <sample>.geneBodyCoverage.txt
│   ├── <sample>.geneBodyCoverage_plot.r
│   ├── <sample>.geneBodyCoverage.pdf
│   ├── <sample>.geneBodyCoverage.curves.r
│   ├── <sample>.geneBodyCoverage.curves.pdf
│   ├── <sample>.geneBodyCoverage.heatmap.r
│   └── <sample>.geneBodyCoverage.heatMap.pdf
├── logs/                                    Per-sample tool logs
│   ├── <sample>.star.log
│   ├── <sample>.samtools.log
│   ├── <sample>.bamcoverage.log
│   └── <sample>.genebody.log
├── .checkpoints/                            Per-step SHA256 checkpoint files
│   ├── <sample>.star.sha256
│   ├── <sample>.deeptools.sha256
│   ├── <sample>.infer.sha256
│   ├── <sample>.rdist.sha256
│   └── <sample>.genebody.sha256
├── annotation.bed12                         Auto-generated BED12 (cached)
├── pipeline_summary.json                    JSON summary of all samples
├── pipeline_summary.tsv                     TSV summary of all samples
└── run_info_latest.json                     Audit trail (config, hashes, versions)
```

---

## STAR Parameters

The following ENCODE-compliant STAR parameters are used:

| Parameter | Value | Purpose |
|-----------|-------|---------|
| `--twopassMode` | `Basic` | 2-pass mapping for novel splice junction discovery |
| `--quantMode` | `TranscriptomeSAM GeneCounts` | Transcriptome BAM + gene-level counts |
| `--outSAMtype` | `BAM SortedByCoordinate` | Coordinate-sorted BAM output |
| `--outSAMstrandField` | `intronMotif` | Strand info for unstranded libraries |
| `--chimSegmentMin` | `15` | Minimum chimeric segment length |
| `--chimJunctionOverhangMin` | `15` | Chimeric junction overhang |
| `--chimScoreMin` | `10` | Minimum chimeric alignment score |
| `--chimScoreDropMax` | `30` | Max score drop for chimeric segments |
| `--chimScoreSeparation` | `10` | Score separation between best chimeric |
| `--chimOutType` | `Junctions SeparateSAMold` | Output chimeric junctions + SAM |
| `--alignSJDBoverhangMin` | `1` | Min overhang for annotated junctions |
| `--alignSJoverhangMin` | `8` | Min overhang for novel junctions |
| `--outFilterMismatchNoverReadLmax` | `0.04` | Max mismatch rate per read length |
| `--alignIntronMin` | `20` | Minimum intron length |
| `--alignIntronMax` | `1000000` | Maximum intron length |
| `--alignMatesGapMax` | `1000000` | Maximum mate pair gap |
| `--limitBAMsortRAM` | auto-detected | RAM for BAM sorting (configurable via `--bam-sort-ram`) |
| `--sjdbGTFfile` | `<GTF>` | Annotation-guided alignment |

---

## Resume and SHA256 Integrity

### How It Works

The pipeline uses **per-step SHA256 checkpoints** for resume awareness. Each
pipeline step writes its own `.sha256` checkpoint file immediately after
successful completion. The checkpoint contains the SHA256 digest of that step's
output files.

On resume, each step's checkpoint is checked:

1. **Checkpoint exists and digest matches current outputs** → step is complete, skip it
2. **Checkpoint missing** → step is incomplete; partial outputs are deleted and the step reruns
3. **Checkpoint exists but digest mismatches** → outputs were corrupted/modified; partial outputs are deleted and the step reruns

This is the **sole resume authority** — no file existence or size checks are used.
Completeness is determined entirely by whether a valid SHA256 checkpoint exists.

### Five Independent Digests

Each sample has up to five checkpoint files:

| Step | Checkpoint File | Files Hashed |
|------|----------------|--------------|
| **STAR** | `<sample>.star.sha256` | Log files, BAI, transcriptome BAM, counts, junctions |
| **deeptools** | `<sample>.deeptools.sha256` | Bigwig file |
| **infer** | `<sample>.infer.sha256` | `strand.txt` |
| **rdist** | `<sample>.rdist.sha256` | `read_distribution.txt` |
| **genebody** | `<sample>.genebody.sha256` | Coverage txt, R scripts, RSeQC PDF, ggplot2 PDFs |

Each file's **name** and **full contents** are fed into a streaming SHA256
hasher. Missing files use a `__MISSING__` sentinel; zero-byte files use a
`__ZERO_BYTES__` sentinel. The digest changes whenever any file is added,
removed, modified, or truncated.

### Checkpoint Format

Each step gets its own file at `.checkpoints/<sample>.<step>.sha256` containing
a single line — the 64-character hex SHA256 digest:

```
a1b2c3d4e5f6789012345678901234567890123456789012345678901234abcd
```

Checkpoints are written atomically (temp file + `fsync` + rename) to prevent
corruption on crashes.

### Resume States

On re-run, the pipeline checks all five steps per sample and determines what
needs to rerun:

| Status | Condition | Action |
|--------|-----------|--------|
| **AllDone** | All 5 checkpoints valid | Skip sample entirely |
| **Phase1Changed** | STAR checkpoint invalid | Clean all outputs, rerun STAR + deeptools + RSeQC |
| **Phase2Changed** | deeptools checkpoint invalid, STAR OK | Clean deeptools + RSeQC outputs, rerun Phase 2 + 3 |
| **Phase3Changed** | Any RSeQC checkpoint invalid, Phases 1-2 OK | Clean only failed RSeQC outputs, rerun Phase 3 |
| **NotDone** | No checkpoints found | Clean partial outputs, process from scratch |

Within Phase 3, each sub-tool (infer, rdist, genebody) has its own checkpoint.
If only genebody fails (e.g., PDF generation error), infer and read_distribution
results are preserved and won't rerun.

### Downstream Invalidation

When a step reruns, all downstream checkpoints are automatically removed:

- Rerunning **STAR** invalidates deeptools + all RSeQC checkpoints
- Rerunning **deeptools** invalidates all RSeQC checkpoints
- Rerunning an **RSeQC tool** only affects that tool's checkpoint

### Resume Example

```bash
# First run — processes all 24 samples
star-rseqc /data/Paired/ -o results

# Second run — skips completed samples
star-rseqc /data/Paired/ -o results
# Output: Resuming: 24/24 samples verified (SHA256 OK), 0 STAR + 0 deeptools + 0 RSeQC to process.

# If genebody PDFs were accidentally deleted for one sample:
star-rseqc /data/Paired/ -o results
# Output: Phase 3 CHANGED: 103N_GBC — cleaning partial QC outputs, will re-run RSeQC only
#         Only genebody reruns; infer + read_distribution checkpoints still valid.
```

---

## Terminal UI

The pipeline features a full-screen alternate-screen terminal interface built
with `crossterm`. Each phase gets its own TUI instance:

```
══════════════════════════════════════════════════════════════════════════════
                             STAR-RSeQC v0.1.0
        STAR 2-Pass Alignment + RSeQC Quality Control | Paired-End RNA-seq
══════════════════════════════════════════════════════════════════════════════
  Phase 1/3 — STAR alignment (24 samples)
  Resumed: 4 sample(s) already completed from previous run
══════════════════════════════════════════════════════════════════════════════
  OVERALL PROGRESS
  ██████████████████████░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░  50%
  12/24 done   Elapsed: 02:34:15   ETA: 02:30:00
  Complete: 20:04:30   Speed: 4.7/min
══════════════════════════════════════════════════════════════════════════════
  ACTIVE JOBS (2/2)
  / 103T_GBC             [STAR alignment] 12m 34s / ~15m 00s
    ████████████████████████████████████████████░░░░░░░░░░░░░░░░░░░░░  83%
  - 104N_CRC             [STAR alignment] 10m 22s / ~15m 00s
    ██████████████████████████████████████░░░░░░░░░░░░░░░░░░░░░░░░░░░  69%
══════════════════════════════════════════════════════════════════════════════
  Completed: 10   Skipped: 4   Failed: 0   Remaining: 10
══════════════════════════════════════════════════════════════════════════════
  RECENT ACTIVITY
    DONE  103N_GBC — STAR alignment
    DONE  103N_GBC — samtools index
    DONE  52N_PACA — STAR alignment
    DONE  52N_PACA — samtools index
    DONE  50T_CRC — STAR alignment
══════════════════════════════════════════════════════════════════════════════
  Ctrl+C to cancel                                       Updated: 17:34:15
```

Features:
- Centered title and subtitle with `═` separator bars
- Phase indicator with resume count from previous runs
- Overall progress bar with percentage, elapsed, ETA, estimated completion time,
  and throughput (samples/min)
- Per-job progress bars with spinning ASCII animation (`|`, `/`, `-`, `\`),
  sample name, current step label, elapsed vs estimated time
- Per-sample bars show indeterminate pulse animation when no average is available
- Color-coded counters: green (completed), yellow (skipped), red (failed),
  white (remaining)
- Color-coded activity log: green (DONE), yellow (SKIP/RESUME), red (FAIL),
  dark red (STOP/cancelled), cyan (INFO)
- Sticky footer with cancel hint and last-updated timestamp
- Graceful `Ctrl+C` — footer changes to red "CANCELLING...", waits for active
  jobs to finish, cleans up partial outputs

---

## Reference Configuration

Default paths are compiled into the binary but can be overridden via command-line
flags or environment variables:

```bash
# Override genome index and GTF
star-rseqc /data/Paired/ \
    --genome-dir /alt/star_index \
    --gtf /alt/annotation.gtf

# Use a different samtools
star-rseqc /data/Paired/ --samtools /usr/bin/samtools

# Use different conda environments
star-rseqc /data/Paired/ \
    --star-env /opt/envs/star \
    --rseqc-env /opt/envs/rseqc \
    --deeptools-env /opt/envs/deeptools

# Environment variables (override compiled defaults, overridden by flags)
export STAR_RSEQC_GENOME_DIR=/alt/star_index
export STAR_RSEQC_GTF=/alt/annotation.gtf
export STAR_RSEQC_STAR_ENV=/opt/envs/star
export STAR_RSEQC_RSEQC_ENV=/opt/envs/rseqc
export STAR_RSEQC_DEEPTOOLS_ENV=/opt/envs/deeptools
```

---

## Examples

```bash
# Basic run on current directory
star-rseqc ./

# Custom output directory and parallelism
star-rseqc /data/transcriptome_01/Paired/ -o batch01-results -j 4 -t 8

# Per-phase job control (more deeptools/RSeQC jobs since they use less RAM)
star-rseqc ./ --star-jobs 2 --deeptools-jobs 6 --rseqc-jobs 6

# QC only on existing BAMs (skip STAR alignment)
star-rseqc ./ --skip-alignment -o existing-results/

# Dry run to check sample discovery and resume status
star-rseqc ./ --dry-run

# Resume after interruption (just re-run the same command)
star-rseqc ./ -o my-results

# Clean up chimeric SAM files after a successful run
star-rseqc ./ --clean-sam
```

### Resource Auto-Detection

Resources are auto-detected from system RAM and CPU count at startup:

- **Jobs**: `(available_RAM - 2 GB) / 38 GB` per STAR job (32 GB genome + 6 GB sort)
- **Threads**: `total_CPUs / jobs`
- **BAM sort RAM**: min(6 GB, 25% usable RAM / jobs)

Override with explicit flags: `-j 2 -t 16 --bam-sort-ram 4000000000`

---

## Architecture

```
main()
 ├── parse_args()              Hand-rolled arg parser (no external crate)
 ├── validate_environment()    Check STAR, samtools, RSeQC, deeptools, Rscript, genome index, GTF
 ├── discover_samples()        Glob *_1P.fastq.gz, pair with *_2P.fastq.gz
 ├── gtf_to_bed12()            Pure-Rust GTF→BED12 converter (cached)
 ├── check_resume_all_parallel()  Per-step SHA256 verification (threaded)
 ├── [dry-run branch]          Print table and exit
 ├── Phase 1: run_work_queue() → run_star_sample()     STAR + samtools index
 │    └── write_step_checkpoint(STEP_STAR)              Checkpoint on success
 ├── Phase 2: run_work_queue() → run_deeptools_phase()  bamCoverage
 │    └── write_step_checkpoint(STEP_DEEPTOOLS)         Checkpoint on success
 └── Phase 3: run_work_queue() → run_rseqc_phase3()    3 tools in parallel
      ├── infer_experiment → write_step_checkpoint(STEP_INFER)
      ├── read_distribution → write_step_checkpoint(STEP_RDIST)
      └── geneBody_coverage2 + ggplot2 PDFs → write_step_checkpoint(STEP_GENEBODY)
```

Key design decisions:
- **No `clap`** — hand-rolled argument parser keeps the dependency tree minimal
- **No `rayon`** — scoped thread work queue with `AtomicUsize` work-stealing
  gives fine-grained control over job slot assignment and TUI updates
- **`crossterm` TUI** — alternate screen with raw mode for a non-blocking
  progress display
- **Streaming SHA256** — files are hashed in 64 KB chunks; multi-GB BAMs are
  excluded from the digest (log files + BAI provide indirect integrity)
- **Per-step checkpoints** — five independent SHA256 digests per sample enable
  granular resume (e.g., only rerun genebody if PDFs were deleted, without
  re-running infer or read_distribution)
- **Atomic checkpoint writes** — temp file + `fsync` + rename prevents
  half-written checkpoints on power loss
- **Downstream invalidation** — rerunning an upstream step automatically removes
  all downstream checkpoints

---

## License

Internal tool. Not yet published under an open-source license.
