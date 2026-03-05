# STAR-RSeQC Installation Architecture

This document describes the three-tier installation system that accommodates different user needs and environments.

## Overview

The STAR-RSeQC distribution package includes three independent installation paths, allowing users to choose based on their needs and existing infrastructure:

```
┌─────────────────────────────────────────────────────────────┐
│                     setup.sh (Main Entry)                   │
└──────────────────┬────────────────────────────────────────────┘
                   │
        ┌──────────┼──────────┐
        │          │          │
        ▼          ▼          ▼
    ┌──────┐  ┌──────┐  ┌──────┐
    │Docker│  │Mamba │  │Manual│
    └──────┘  └──────┘  └──────┘
```

---

## Path 1: Docker Installation (Fastest)

**Best for:** Users who want everything pre-configured with minimal setup.

### Flow

```
1. Check if Docker is installed
   ├─ YES → Build/run Docker image
   └─ NO  → Ask user: Auto-install or Manual?
           ├─ Auto → Install Docker (ubuntu/debian/fedora/rhel)
           │         └─ Build image, create docker-run.sh
           └─ Manual → Provide links to Docker docs
                       └─ User manually installs Docker
```

### Files Used

- `Dockerfile` - Multi-stage build (Rust builder + mamba runtime)
- `docker-compose.yml` - Docker Compose configuration
- `.dockerignore` - Build optimization

### Features

- ✅ **Zero dependency installation** — everything in one container
- ✅ **Reproducible environment** — same setup everywhere
- ✅ **Portable** — runs on any system with Docker
- ✅ **Auto Docker install** — script can install Docker for you (Ubuntu, Debian, Fedora, RHEL)
- ⏱️ **Time**: 5-15 minutes (depending on Docker availability)

### Advantages

- All dependencies pre-compiled
- No conda/mamba needed
- Works on macOS, Linux, Windows (with WSL)
- Easy to remove (`docker rmi star-rseqc:latest`)

### Usage

```bash
./setup.sh
# Select: 1 (Docker)
# [optional: auto-install Docker]
# ./docker-run.sh /path/to/fastq
```

---

## Path 2: Automatic Mamba Installation

**Best for:** Users who prefer conda/mamba and want fully automated setup.

### Flow

```
1. Check for existing conda/mamba
   ├─ YES → Create STAR and RSeQC environments
   └─ NO  → Ask user: Auto-install Miniforge or Manual?
           ├─ Auto → Download & install Miniforge
           │         └─ Create environments
           │         └─ Build binary
           │         └─ Write config.json
           └─ Manual → Provide conda installation links
```

### Installation Steps

1. **Auto-detect or install conda**
   - Searches: `~/miniforge3`, `~/mambaforge`, `~/miniconda3`, `/opt/conda`
   - If not found: downloads Miniforge installer (~150 MB)

2. **Create mamba environments**
   ```bash
   mamba create -n star -c bioconda -c conda-forge star=2.7.11b samtools
   mamba create -n rseqc -c bioconda -c conda-forge rseqc=5.0.4 python
   ```

3. **Build binary**
   ```bash
   cargo build --release
   cp target/release/star-rseqc ~/.local/bin/
   ```

4. **Write config**
   ```json
   {
     "genome_dir": "/user/provided/path",
     "gtf": "/user/provided/path",
     "star_env": "$HOME/miniforge3/envs/star",
     "rseqc_env": "$HOME/miniforge3/envs/rseqc",
     "samtools": "$HOME/miniforge3/envs/star/bin/samtools"
   }
   ```

### Features

- ✅ **Auto-installs Miniforge** if conda not present
- ✅ **Creates isolated environments** (STAR + RSeQC separate)
- ✅ **Fully automated** — no manual steps needed
- ✅ **Config file generated** automatically
- ⏱️ **Time**: 20-30 minutes (includes 2-3 GB download)

### Advantages

- Uses conda/mamba (familiar to bioinformaticians)
- Can leverage existing conda packages
- Conda environments easy to recreate or modify
- Better for HPC systems with module support

### Usage

```bash
./setup.sh
# Select: 2 (Automatic)
# Follow prompts (mostly just press Enter)
# Enter reference file paths when asked
```

---

## Path 3: Manual Installation

**Best for:** Advanced users, HPC clusters, custom environments.

### Provides

- Links to official installation pages:
  - [Miniforge](https://github.com/conda-forge/miniforge)
  - [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
  - [Docker](https://docs.docker.com/engine/install/)

- Step-by-step instructions for:
  1. Installing conda/mamba
  2. Creating STAR environment
  3. Creating RSeQC environment
  4. Building STAR-RSeQC from source
  5. Creating config.json manually

### Features

- ✅ **Full control** over each step
- ✅ **Flexible** — use system packages if available
- ✅ **HPC-friendly** — can use module-loaded tools
- ✅ **Debugging** — easier to troubleshoot custom setups
- ⏱️ **Time**: 30-60 minutes (depends on user experience)

### Advantages

- Maximum control over configuration
- Can use system-wide STAR/RSeQC installations
- Integrates with existing conda environments
- Better for shared/HPC systems

### Usage

```bash
./setup.sh
# Select: 3 (Manual)
# Read provided links
# Follow manual steps
# Create ~/.config/star-rseqc/config.json manually
```

---

## Configuration Resolution

All three paths produce the same end result: a working binary with proper configuration.

```
Priority Order (highest to lowest):
1. CLI flags               (--genome-dir, --star-env, etc.)
2. ~/.config/star-rseqc/config.json
3. Auto-detection          (find_conda_env, find_samtools)
4. Empty PathBuf           (caught by validation)
```

---

## Decision Tree

```
Choose installation method based on:

Does Docker appeal to you?
├─ YES: Use Path 1 (Docker)
│       Fast, zero config, fully containerized
│
├─ NO: Do you want automatic installation?
│      ├─ YES: Use Path 2 (Automatic Mamba)
│      │        ~25 minutes, installs everything for you
│      │
│      └─ NO: Use Path 3 (Manual)
│             ~45 minutes, full control, step-by-step
│
Do you have >20 GB internet bandwidth?
├─ YES: Any path works well
└─ NO: Docker is fastest (~500 MB)
       Manual might be better (reuse existing installs)
```

---

## Package Distribution

The `create-package.sh` script bundles everything:

```bash
./create-package.sh [VERSION]
# Creates: STAR-RSeQC-0.2.0.zip

# Contents:
#   ├── star-rseqc (binary)
#   ├── setup.sh (main installer)
#   ├── docker/ (Dockerfile, etc.)
#   ├── docs/ (README, Quick Start, Troubleshooting)
#   ├── INSTALL.txt (installation guide)
#   ├── VERSION.txt (version info)
#   └── README_PACKAGE.txt (package contents)
```

**Distribution:**
1. Users download `STAR-RSeQC-0.2.0.zip`
2. Extract: `unzip STAR-RSeQC-0.2.0.zip`
3. Run: `bash setup.sh`
4. Choose their preferred path

---

## Comparison Table

| Aspect | Docker | Auto Mamba | Manual |
|--------|--------|-----------|--------|
| **Setup Time** | 5-15 min | 20-30 min | 30-60 min |
| **Dependencies** | Docker only | Bash + Internet | Bash + Internet + Build tools |
| **Internet Required** | 500 MB | 2-3 GB | 2-3 GB (or reuse existing) |
| **Portability** | Excellent | Good | Good |
| **Customization** | Limited | Moderate | Full |
| **HPC-Friendly** | No | Yes | Yes |
| **Requires Root?** | Docker install needs sudo | No | Maybe (mamba install) |
| **Works Offline?** | No | No | Yes (if tools already installed) |

---

## Troubleshooting by Path

### Docker Path

**Issue**: "Docker daemon not running"
- Solution: `sudo systemctl start docker`

**Issue**: "Permission denied while trying to connect"
- Solution: `sudo usermod -aG docker $USER && newgrp docker`

### Mamba Path

**Issue**: "conda not found after installation"
- Solution: `source ~/.bashrc && conda --version`

**Issue**: "STAR environment not created"
- Solution: `mamba env list` to verify, re-run setup if missing

### Manual Path

**Issue**: "STAR not found in expected path"
- Solution: Manually specify with `--star-env /custom/path`

**Issue**: "samtools command not found"
- Solution: Ensure it's in PATH or specify with `--samtools /full/path`

---

## Advanced: Custom Docker Image

Build your own image with custom reference files:

```dockerfile
FROM star-rseqc:latest

# Copy reference files into image
COPY hg38_star_index /ref/genome
COPY hg38.gtf /ref/annotation.gtf

# Update config
RUN cat > /root/.config/star-rseqc/config.json << EOF
{
  "genome_dir": "/ref/genome",
  "gtf": "/ref/annotation.gtf",
  "star_env": "/opt/conda/envs/star",
  "rseqc_env": "/opt/conda/envs/rseqc",
  "samtools": "/opt/conda/envs/star/bin/samtools"
}
EOF
```

Build: `docker build -t star-rseqc:custom .`

---

## Summary

- **Docker**: Fastest, most isolated (if Docker available)
- **Auto Mamba**: Balanced, fully automatic
- **Manual**: Maximum control, HPC-friendly

All three produce identical functionality. Choose based on your environment and comfort level.
