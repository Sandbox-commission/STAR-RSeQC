# star-rseqc

**STAR 2-pass alignment + RSeQC quality control pipeline for paired-end RNA-seq**

A high-performance, resume-aware pipeline written in Rust that automates STAR
two-pass alignment, BAM indexing, and RSeQC quality control for bulk paired-end
RNA-seq experiments. Features a full-screen terminal UI with real-time progress
tracking and directory-wise SHA256 integrity verification.

---

## Quick Install

```bash
curl -fsSL https://raw.githubusercontent.com/Sandbox-commission/STAR-RSeQC/main/setup.sh | bash
```

This interactive script will:
- Auto-detect your conda/mamba installation
- Create `star` and `rseqc` conda environments
- Prompt for reference file paths (or build a STAR index if needed)
- Install the pre-built binary (or compile from source)
- Write configuration to `~/.config/star-rseqc/config.json`

After setup, add `~/.local/bin` to your PATH and run:
```bash
star-rseqc /path/to/fastq/directory
```

---

## Table of Contents

- [Quick Install](#quick-install)
- [Features](#features)
- [Requirements](#requirements)
- [Getting Reference Files](#getting-reference-files)
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
- [Troubleshooting](#troubleshooting)
- [Example Output](#example-output)
  - [Understanding your QC results](#understanding-your-qc-results)
- [Architecture](#architecture)
- [References](#references)
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
  activity log (built with `crossterm`)
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

| Tool | Version | Setup |
|------|---------|--------|
| [STAR](https://github.com/alexdobin/STAR) | v2.7.11b+ | Auto-installed via `setup.sh` in `star` conda environment |
| [samtools](http://www.htslib.org/) | v1.15+ | Auto-installed via `setup.sh` in `star` conda environment |
| [RSeQC](http://rseqc.sourceforge.net/) | v5.0+ | Auto-installed via `setup.sh` in `rseqc` conda environment |
| conda/mamba | — | Required; auto-detected from standard install locations |

### Reference files

| File | Setup |
|------|--------|
| STAR genome index | Specified interactively during `setup.sh`; or auto-built from FASTA+GTF if needed |
| GTF annotation | Specified interactively during `setup.sh` |
| Reference FASTA | Required only if building STAR index from scratch (optional) |

All paths are stored in `~/.config/star-rseqc/config.json` and auto-detected or manually overridable via CLI flags.

---

## Getting Reference Files

If you're new to bioinformatics, you'll need two things before running the pipeline:
a **GTF annotation file** (tells STAR where genes are) and a **STAR genome index**
(a pre-processed version of the genome that STAR can align reads against).

### Step 1: Download genome and annotation

Pick the organism you're working with. For **human (GRCh38)**:

```bash
# Create a directory for reference files
mkdir -p ~/references && cd ~/references

# Download the genome FASTA (DNA sequence, ~900 MB compressed)
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Download the GTF annotation (gene locations, ~50 MB compressed)
wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
gunzip Homo_sapiens.GRCh38.113.gtf.gz
```

For **mouse (GRCm39)**, replace `homo_sapiens` with `mus_musculus` and
`GRCh38` with `GRCm39` in the URLs above. For other organisms, browse
https://ftp.ensembl.org/pub/release-113/.

### Step 2: Build the STAR genome index

This step takes 30-60 minutes and needs ~32 GB of RAM for the human genome.
You only need to do this once.

```bash
# Activate the STAR conda environment
conda activate star

# Build the index (adjust --sjdbOverhang to your read length minus 1;
# 100 works for most Illumina experiments with 101 bp reads)
mkdir -p ~/references/star_index
STAR --runMode genomeGenerate \
     --genomeDir ~/references/star_index \
     --genomeFastaFiles ~/references/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
     --sjdbGTFfile ~/references/Homo_sapiens.GRCh38.113.gtf \
     --sjdbOverhang 100 \
     --runThreadN 8

conda deactivate
```

### Step 3: Run the pipeline

```bash
star-rseqc /path/to/your/fastq/files \
    --genome-dir ~/references/star_index \
    --gtf ~/references/Homo_sapiens.GRCh38.113.gtf
```

Or save these paths permanently so you don't need to type them every time:

```bash
mkdir -p ~/.config/star-rseqc
cat > ~/.config/star-rseqc/config.json << EOF
{
  "genome_dir": "$HOME/references/star_index",
  "gtf": "$HOME/references/Homo_sapiens.GRCh38.113.gtf"
}
EOF

# Now just run:
star-rseqc /path/to/your/fastq/files
```

---

## Installation

### Automated Setup (Recommended)

The easiest way to get started is to use the provided setup script (see [Quick Install](#quick-install)):

```bash
curl -fsSL https://raw.githubusercontent.com/Sandbox-commission/STAR-RSeQC/main/setup.sh | bash
```

### Manual Install

1. **From source:**
   ```bash
   git clone https://github.com/Sandbox-commission/STAR-RSeQC.git
   cd STAR-RSeQC
   cargo build --release
   cp target/release/star-rseqc ~/.local/bin/
   ```

2. **Create conda environments:**
   ```bash
   mamba create -n star -c bioconda -c conda-forge star=2.7.11b samtools
   mamba create -n rseqc -c bioconda -c conda-forge rseqc=5.0.4 python
   ```

3. **Create config file:**
   ```bash
   mkdir -p ~/.config/star-rseqc
   cat > ~/.config/star-rseqc/config.json << 'EOF'
   {
     "genome_dir": "/path/to/star/index",
     "gtf": "/path/to/annotation.gtf",
     "star_env": "/path/to/miniforge3/envs/star",
     "rseqc_env": "/path/to/miniforge3/envs/rseqc",
     "samtools": "/path/to/samtools"
   }
   EOF
   ```

4. **Add to PATH:**
   ```bash
   export PATH="$HOME/.local/bin:$PATH"
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
SAMPLE1_CONTROL_1P.fastq.gz   →  sample = SAMPLE1_CONTROL
SAMPLE1_CASE_1P.fastq.gz    →  sample = SAMPLE1_CASE
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
# Output: RSeQC changed: SAMPLE1_CONTROL (was fedcba987654..., now 1234abcd5678...)
#         1 sample(s) have corrupted/changed outputs — will re-process.
```

---

## Terminal UI

The pipeline features a full-screen terminal interface built with `crossterm`:

```
╔══════════════════════════════════════════════════════════════════╗
║  STAR-RSeQC v0.2.0                                               ║
║  STAR 2-Pass Alignment + RSeQC Quality Control | Paired-End      ║
╠══════════════════════════════════════════════════════════════════╣
║  STAR 2-pass + RSeQC (24 samples, 2x16t [auto])                 ║
║                                                                  ║
║  OVERALL PROGRESS                                                ║
║  ██████████████░░░░░░░░░░░░░░░░  50%                             ║
║  12/24 done   Elapsed: 02:34:15   ETA: 02:30:00   Speed: 0.1/min║
╠══════════════════════════════════════════════════════════════════╣
║  ACTIVE JOBS (2/2)                                               ║
║  | SAMPLE1_CONTROL       [STAR alignment] 12:34 / ~25:00        ║
║    ████████████████░░░░░░░░░░░░░░░░░░░░  50%                     ║
║  / SAMPLE1_CASE          [RSeQC: gene body coverage] 02:11      ║
║    ████████░░░░░░░░░░░░░░░░░░░░░░░░░░░░   8%                     ║
╠══════════════════════════════════════════════════════════════════╣
║  Completed: 10   Skipped: 2   Failed: 0   Remaining: 12         ║
╠══════════════════════════════════════════════════════════════════╣
║  RECENT ACTIVITY                                                 ║
║  DONE  SAMPLE1_CONTROL — star:a1b2c3d4e5f6… rseqc:7890abcdef12… ║
║  DONE  SAMPLE1_CONTROL — STAR alignment                         ║
║  DONE  SAMPLE1_CONTROL — samtools index                         ║
║                                                                  ║
║  Ctrl+C to cancel                            Updated: 14:32:15  ║
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

Configuration values are resolved in a 3-layer priority system:

1. **CLI flags** (highest priority)
   - `--genome-dir`, `--gtf`, `--star-env`, `--rseqc-env`, `--samtools`
   - Special: `--star-env auto` or `--rseqc-env auto` trigger automatic environment search

2. **Config file** (`~/.config/star-rseqc/config.json`)
   - User-local configuration persisted by `setup.sh`
   - Example:
     ```json
     {
       "genome_dir": "/mnt/reference/star_hg38_101bp",
       "gtf": "/mnt/reference/Homo_sapiens.GRCh38.113.gtf",
       "star_env": "/home/user/miniforge3/envs/star",
       "rseqc_env": "/home/user/miniforge3/envs/rseqc",
       "samtools": "/home/user/miniforge3/envs/star/bin/samtools"
     }
     ```

3. **Auto-detection** (lowest priority)
   - Searches for `star` and `rseqc` conda environments in standard locations:
     - `~/miniforge3/envs/{name}`
     - `~/mambaforge/envs/{name}`
     - `~/miniconda3/envs/{name}`
     - `/opt/conda/envs/{name}`
   - Searches for `samtools` in STAR environment `bin/` and system PATH

### Usage Examples

```bash
# Use config file + auto-detection (no flags needed)
star-rseqc /data/Paired/

# Override with CLI flags (highest priority)
star-rseqc /data/Paired/ \
    --genome-dir /alt/star_index \
    --gtf /alt/annotation.gtf

# Auto-detect conda environments
star-rseqc /data/Paired/ \
    --star-env auto \
    --rseqc-env auto

# Mix config file, auto-detection, and CLI overrides
star-rseqc /data/Paired/ \
    --samtools /usr/local/bin/samtools \
    --threads 16
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

## Troubleshooting

### STAR not found

```
STAR binary not found: /home/user/miniforge3/envs/star/bin/STAR
```

The `star` conda environment is missing or `--star-env` points to the wrong location. Fix:
```bash
mamba create -n star -c bioconda -c conda-forge star=2.7.11b samtools
star-rseqc ./ --star-env auto
```

### samtools version too old

```
samtools version 1.9 detected; version 1.15+ is recommended.
```

Upgrade samtools in the STAR environment:
```bash
mamba install -n star -c bioconda -c conda-forge "samtools>=1.15"
```

### Docker socket errors

```
Cannot connect to the Docker daemon
```

If using the system Docker daemon (not Docker Desktop):
```bash
export DOCKER_HOST="unix:///var/run/docker.sock"
```
Or ensure Docker is running: `sudo systemctl start docker`

### Resume behavior

The pipeline uses SHA256 checksums over output files to detect completed samples. If outputs are accidentally deleted or corrupted, the pipeline will detect the change and re-process those samples automatically. To force a full re-run, delete the `.checkpoints/` directory inside the output folder.

### GTF conversion failures

```
GTF->BED12 produced zero transcripts
```

This usually means the GTF file is empty, truncated, or uses a non-standard format. Verify your GTF file contains lines with feature type `exon` and a `transcript_id` attribute. Alternatively, provide a pre-built BED12 file with `--bed`.

### Permission errors

If output files can't be written, check that you have write permissions on the output directory. In Docker, ensure the mounted volumes have correct permissions and the `user:` setting in `docker-compose.yml` matches your host UID/GID (default `1000:1000`).

---

## Example Output

### pipeline_summary.json

After a successful run, `pipeline_summary.json` contains one entry per sample:

```json
[
  {
    "sample": "SAMPLE1_CONTROL",
    "sha256": "star:a1b2c3d4e5f67890|rseqc:fedcba9876543210",
    "log_final": true,
    "bam_sorted": true,
    "bam_index": true,
    "bam_transcriptome": true,
    "gene_counts": true,
    "splice_junctions": true,
    "chimeric_junction": true,
    "chimeric_sam": true,
    "strand_qc": true,
    "genebody_txt": true,
    "genebody_r": true,
    "genebody_curves_pdf": true,
    "genebody_heatmap_pdf": true,
    "readdist_qc": true
  }
]
```

### Understanding your QC results

After the pipeline finishes, check the files in the `qc/` folder. Here's what
each one tells you and what "good" vs "bad" looks like:

#### Strandedness (`*.strand.txt`)

This tells you whether your RNA library was prepared with strand information.

| Result | Meaning |
|--------|---------|
| "1++,1--,2+-,2-+" fraction **> 0.8** | **Stranded** (sense) — this is normal for most modern kits (e.g. Illumina TruSeq Stranded) |
| "1+-,1-+,2++,2--" fraction **> 0.8** | **Reverse-stranded** — also normal, depends on the kit used |
| Both fractions **near 0.5** | **Unstranded** — older protocols or some poly-A kits; still fine for gene-level analysis |

If you're unsure what to expect, ask whoever prepared the RNA library which kit they used.

#### Gene body coverage (`*.geneBodyCoverage.txt` and `*.curves.pdf`)

Open the PDF — it shows a curve from the 5' end (start) to the 3' end (end) of genes.

| Shape | Meaning |
|-------|---------|
| **Flat/even curve** | Good — RNA was intact when sequenced |
| **Strong drop-off at the 5' end** (left side lower) | RNA may be degraded — the 3' end is preserved but the 5' end is lost. Common in old or poorly stored samples |
| **Spike at one end only** | Possible bias in library preparation |

Moderate 3' bias is common and usually acceptable. A severe drop (5' signal less than half of 3') suggests degradation that may affect downstream analysis.

#### Read distribution (`*.read_distribution.txt`)

This shows where your reads mapped across the genome.

| Region | Typical good result | Concern if... |
|--------|-------------------|---------------|
| **CDS exons** | 40-60% of reads | Very low (< 20%) — possible DNA contamination |
| **UTR exons (5' and 3')** | 15-30% of reads | — |
| **Introns** | 10-25% of reads | Very high (> 50%) — may indicate DNA contamination or pre-mRNA |
| **Intergenic** | < 10% of reads | Very high (> 30%) — possible DNA contamination or poor rRNA depletion |

**In plain language:** most reads should land on known gene regions (exons).
If a large fraction maps between genes (intergenic) or within gene bodies but
not on exons (intronic), something may be off with the sample or library prep.

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
- **No `clap`** — hand-rolled argument parser for minimal dependencies
- **No `rayon`** — scoped thread work queue with `AtomicUsize` work-stealing
  gives fine-grained control over job slot assignment and TUI updates
- **`crossterm` TUI** — alternate screen with raw mode for clean terminal rendering
- **Streaming SHA256** — files are hashed in 64 KB chunks to handle large BAMs
  without loading them into memory
- **Directory-wise digests** — two SHA256 hashes per sample (STAR + RSeQC)
  instead of one combined hash, enabling precise identification of which output
  directory was corrupted

---

## References

If you use STAR-RSeQC in your research, please cite the underlying tools:

- **STAR**: Dobin A, Davis CA, Schlesinger F, et al. *STAR: ultrafast universal RNA-seq aligner.* Bioinformatics. 2013;29(1):15-21. doi:[10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635) | [PMID: 23104886](https://pubmed.ncbi.nlm.nih.gov/23104886/)

- **RSeQC**: Wang L, Wang S, Li W. *RSeQC: quality control of RNA-seq experiments.* Bioinformatics. 2012;28(16):2184-2185. doi:[10.1093/bioinformatics/bts356](https://doi.org/10.1093/bioinformatics/bts356) | [PMID: 22743226](https://pubmed.ncbi.nlm.nih.gov/22743226/)

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
