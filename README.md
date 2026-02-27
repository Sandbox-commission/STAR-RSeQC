# star-rseqc

**STAR 2-pass alignment + RSeQC quality control pipeline for paired-end RNA-seq**

A high-performance, resume-aware pipeline written in Rust that automates STAR
two-pass alignment, BAM indexing, and RSeQC quality control for bulk paired-end
RNA-seq experiments. Features a full-screen terminal UI with real-time progress
tracking and directory-wise SHA256 integrity verification.

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
- [Pipeline Steps](#pipeline-steps)
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

- **STAR 2-pass alignment** with chimeric junction detection, transcriptome BAM,
  and gene-level quantification (ENCODE-compliant parameters)
- **RSeQC quality control**: strandedness inference, gene body coverage, and
  read distribution analysis
- **Pure-Rust GTF-to-BED12 conversion** — no external tools like `gtfToGenePred`
  needed; the annotation is converted automatically on first run and cached
- **Full-screen TUI** — real-time progress monitor with per-sample spinners,
  overall progress bar, active job slots, elapsed/ETA timers, and a scrolling
  activity log (built with `crossterm`, matching the style of
  [msi-calc](https://github.com/your-repo/msi-calc))
- **Directory-wise SHA256 resume** — on completion, separate SHA256 digests are
  computed for STAR outputs (10 files) and RSeQC outputs (6 files); on re-run,
  digests are verified and only corrupted or incomplete samples are re-processed
- **Parallel execution** — configurable number of concurrent sample jobs, each
  with its own thread allocation for STAR
- **Graceful cancellation** — `Ctrl+C` signals all running jobs to stop cleanly;
  partial STAR outputs are removed so corrupted BAMs never persist
- **Dry-run mode** — list discovered samples and their resume status without
  executing anything

---

## Requirements

| Tool | Version | Conda Environment |
|------|---------|-------------------|
| [STAR](https://github.com/alexdobin/STAR) | v2.7.11b+ | `/home/cml/miniforge3/envs/star` |
| [samtools](http://www.htslib.org/) | v1.15+ | `/home/cml/Downloads/samtools-1.15.1/samtools` |
| [RSeQC](http://rseqc.sourceforge.net/) | v5.0+ | `/home/cml/miniforge3/envs/RSeQC` |
| Rust toolchain | 1.70+ (edition 2021) | — |

### Reference files

| File | Default Path |
|------|-------------|
| STAR genome index | `/home/cml/humandb/transcriptomeindex/ensembl113/star_hg38_101bp_index` |
| GTF annotation | `/home/cml/humandb/transcriptomeindex/ensembl113/Homo_sapiens.GRCh38.113.gtf` |
| Reference FASTA | `/home/cml/humandb/transcriptomeindex/ensembl113/genome.fa` |

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
| `-o, --output <DIR>` | Output directory for all results | `star-rseqc-results` |
| `-j, --jobs <N>` | Number of samples processed in parallel | `2` |
| `-t, --threads <N>` | Threads allocated per STAR alignment job | `16` |
| `--genome-dir <DIR>` | STAR genome index directory | (see defaults) |
| `--gtf <FILE>` | GTF annotation file | (see defaults) |
| `--bed <FILE>` | Pre-built BED12 file for RSeQC (auto-generated from GTF if omitted) | auto |
| `--samtools <PATH>` | Path to samtools binary | (see defaults) |
| `--star-env <DIR>` | STAR conda environment prefix | (see defaults) |
| `--rseqc-env <DIR>` | RSeQC conda environment prefix | (see defaults) |
| `--bam-sort-ram <BYTES>` | RAM limit for STAR BAM sorting | `30000000000` (30 GB) |
| `--skip-qc` | Skip all RSeQC steps (alignment only) | off |
| `--skip-alignment` | Skip STAR alignment (run QC on existing BAMs) | off |
| `--dry-run` | List discovered samples and resume status without executing | off |
| `-h, --help` | Print the full help message | — |

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

Both R1 and R2 must exist for a sample to be included. Samples without a
matching R2 are skipped with a warning.

---

## Pipeline Steps

For each sample, the pipeline executes three stages sequentially:

### Step 1: STAR 2-Pass Alignment

Runs STAR in two-pass mode with chimeric junction detection:

- Produces coordinate-sorted BAM (`*_Aligned.sortedByCoord.out.bam`)
- Produces transcriptome BAM (`*_Aligned.toTranscriptome.out.bam`)
- Generates gene-level counts (`*_ReadsPerGene.out.tab`)
- Detects chimeric reads for fusion discovery (`*_Chimeric.out.junction`)
- Logs STAR stdout/stderr to `logs/<sample>.star.log`

If STAR fails, all partial output files for that sample (including any `_STARtmp`
directory) are cleaned up automatically.

### Step 2: samtools index

Indexes the coordinate-sorted BAM to produce `*.bam.bai`, required by all
downstream tools.

### Step 3: RSeQC Quality Control

Three RSeQC modules are run (each is non-fatal — failure of one does not block
the others):

| Module | Output | Purpose |
|--------|--------|---------|
| `infer_experiment.py` | `<sample>.strand.txt` | Library strandedness (sense/antisense/unstranded) |
| `geneBody_coverage.py` | `<sample>.geneBodyCoverage.{txt,r,curves.pdf,heatMap.pdf}` | 5'-to-3' coverage uniformity |
| `read_distribution.py` | `<sample>.read_distribution.txt` | Read distribution across genomic features |

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
├── qc/                                      RSeQC quality control outputs
│   ├── <sample>.strand.txt
│   ├── <sample>.geneBodyCoverage.txt
│   ├── <sample>.geneBodyCoverage.r
│   ├── <sample>.geneBodyCoverage.curves.pdf
│   ├── <sample>.geneBodyCoverage.heatMap.pdf
│   └── <sample>.read_distribution.txt
├── logs/                                    Per-sample STAR stderr logs
│   └── <sample>.star.log
├── .checkpoints/                            SHA256 checkpoint files
│   └── <sample>.sha256
├── annotation.bed12                         Auto-generated BED12 (cached)
├── pipeline_summary.json                    JSON summary of all samples
└── pipeline_summary.tsv                     TSV summary of all samples
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
| `--limitBAMsortRAM` | `30000000000` | 30 GB RAM for BAM sorting (configurable) |
| `--sjdbGTFfile` | `<GTF>` | Annotation-guided alignment |

---

## Resume and SHA256 Integrity

### How It Works

The pipeline uses **directory-wise SHA256 digests** to verify output integrity
and enable safe resume. Instead of hashing multi-GB input FASTQ files (which
would be prohibitively slow), it computes digests over the **output files** after
each sample completes.

Two independent SHA256 digests are computed per sample:

1. **STAR digest** — covers 10 files in the `star/` directory:
   - `*_Log.out`, `*_Log.progress.out`, `*_Log.final.out`
   - `*_Aligned.sortedByCoord.out.bam`, `*_Aligned.sortedByCoord.out.bam.bai`
   - `*_Aligned.toTranscriptome.out.bam`
   - `*_ReadsPerGene.out.tab`, `*_SJ.out.tab`
   - `*_Chimeric.out.junction`, `*_Chimeric.out.sam`

2. **RSeQC digest** — covers 6 files in the `qc/` directory:
   - `*.strand.txt`
   - `*.geneBodyCoverage.txt`, `*.geneBodyCoverage.r`
   - `*.geneBodyCoverage.curves.pdf`, `*.geneBodyCoverage.heatMap.pdf`
   - `*.read_distribution.txt`

Each file's **name** and **full contents** are fed into a streaming SHA256
hasher. If a file is missing, a `__MISSING__` sentinel is hashed instead, so the
digest changes whenever a file is added, removed, or modified.

### Checkpoint Format

Each sample gets exactly **one** checkpoint file at
`.checkpoints/<sample>.sha256`:

```
star:a1b2c3d4e5f6789012345678901234567890123456789012345678901234abcd
rseqc:fedcba9876543210fedcba9876543210fedcba9876543210fedcba98765432ef
```

Two lines, each prefixed with `star:` or `rseqc:` followed by the 64-character
hex SHA256 digest.

### Resume States

On re-run, the pipeline recomputes digests from files on disk and compares
against the stored checkpoint:

| Status | Condition | Action |
|--------|-----------|--------|
| **DONE (SHA256 OK)** | Both STAR and RSeQC digests match | Skip sample |
| **STAR CHANGED** | Only the STAR digest differs | Re-process entire sample |
| **RSeQC CHANGED** | Only the RSeQC digest differs | Re-process entire sample |
| **STAR+RSeQC CHANGED** | Both digests differ | Re-process entire sample |
| **PENDING** | No checkpoint file exists | Process sample |

This design provides clear diagnostics: you can immediately tell whether
alignment outputs, QC outputs, or both were corrupted or deleted.

### Resume Example

```bash
# First run — processes all 24 samples
star-rseqc /data/Paired/ -o results

# Second run — skips completed samples (SHA256 verified)
star-rseqc /data/Paired/ -o results
# Output: Resuming: 24/24 samples verified (SHA256 OK), 0 to process.

# If some QC files were accidentally deleted:
star-rseqc /data/Paired/ -o results
# Output: RSeQC changed: 103N_GBC (was fedcba987654..., now 1234abcd5678...)
#         1 sample(s) have corrupted/changed outputs — will re-process.
```

---

## Terminal UI

The pipeline features a full-screen terminal interface built with `crossterm`:

```
╔══════════════════════════════════════════════════════════════════╗
║  star-rseqc v0.1.0  │  STAR + RSeQC RNA-seq Pipeline           ║
╠══════════════════════════════════════════════════════════════════╣
║  Phase: Processing samples                                      ║
║                                                                  ║
║  Overall: ██████████████░░░░░░░░░░░░░░░░  12/24  (50.0%)       ║
║  Elapsed: 02:34:15  │  ETA: ~02:30:00                           ║
║                                                                  ║
║  ── Active Jobs ──                                               ║
║  ⠋ [1] 103T_GBC     STAR alignment          (12m 34s)          ║
║  ⠹ [2] 104N_CRC     RSeQC: gene body cov.   (02m 11s)         ║
║                                                                  ║
║  Completed: 10  │  Skipped: 2  │  Failed: 0                    ║
║                                                                  ║
║  ── Recent Activity ──                                           ║
║    DONE  103N_GBC — star:a1b2c3d4e5f6… rseqc:7890abcdef12…     ║
║    DONE  52N_PACA — STAR alignment                               ║
║    DONE  52N_PACA — samtools index                               ║
║                                                                  ║
║  Press Ctrl+C to cancel gracefully                               ║
╚══════════════════════════════════════════════════════════════════╝
```

Features:
- Unicode box-drawing characters with color-coded status
- Animated spinners (braille pattern) for active jobs
- Per-job step labels showing exactly what each slot is doing
- Overall progress bar with percentage and ETA
- Scrolling activity log showing recent completions/failures
- Graceful `Ctrl+C` handling — waits for active jobs to finish

---

## Reference Configuration

Default paths are compiled into the binary. Override any of them with command-line
flags:

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
    --rseqc-env /opt/envs/rseqc
```

---

## Examples

```bash
# Basic run on current directory
star-rseqc ./

# Custom output directory and parallelism
star-rseqc /data/transcriptome_01/Paired/ -o batch01-results -j 4 -t 8

# Alignment only (skip all QC)
star-rseqc ./ --skip-qc

# QC only on existing BAMs (skip STAR alignment)
star-rseqc ./ --skip-alignment -o existing-results/

# Dry run to check sample discovery and resume status
star-rseqc ./ --dry-run

# Resume after interruption (just re-run the same command)
star-rseqc ./ -o my-results
```

### Resource Planning

Total CPU usage = `jobs` x `threads`:

| Jobs | Threads | Total Cores | RAM (approx) |
|------|---------|-------------|--------------|
| 2 | 16 | 32 | ~64 GB |
| 4 | 8 | 32 | ~128 GB |
| 4 | 16 | 64 | ~128 GB |

Each STAR job loads the full genome index (~32 GB for human) into shared memory.
Multiple jobs share this memory, but each needs `--bam-sort-ram` (default 30 GB)
for BAM sorting. Adjust `-j` and `-t` based on your system's available cores and
RAM.

---

## Architecture

```
main()
 ├── parse_args()              Hand-rolled arg parser (no external crate)
 ├── validate_environment()    Check STAR, samtools, RSeQC, genome index, GTF
 ├── discover_samples()        Glob *_1P.fastq.gz, pair with *_2P.fastq.gz
 ├── gtf_to_bed12()            Pure-Rust GTF→BED12 converter (cached)
 ├── check_resume() per sample SHA256 output verification
 ├── [dry-run branch]          Print table and exit
 └── run_work_queue()          Scoped thread pool with atomic work-stealing
      ├── DisplayThread        Separate rendering thread (crossterm, 100ms tick)
      └── process_sample()     Per-sample pipeline:
           ├── STAR 2-pass     run_cancellable() with cleanup on failure
           ├── samtools index   run_cancellable()
           └── RSeQC (3 tools) infer_experiment, geneBody_coverage, read_distribution
                                Non-fatal: failures logged but don't abort sample
```

Key design decisions:
- **No `clap`** — hand-rolled argument parser matches the style of `msi-calc`
- **No `rayon`** — scoped thread work queue with `AtomicUsize` work-stealing
  gives fine-grained control over job slot assignment and TUI updates
- **`crossterm` TUI** — alternate screen with raw mode, matching `msi-calc`'s
  rendering pattern
- **Streaming SHA256** — files are hashed in 64 KB chunks to handle large BAMs
  without loading them into memory
- **Directory-wise digests** — two SHA256 hashes per sample (STAR + RSeQC)
  instead of one combined hash, enabling precise identification of which output
  directory was corrupted

---

## License

Internal tool. Not yet published under an open-source license.
