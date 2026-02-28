#!/bin/bash
set -euo pipefail

# STAR-RSeQC Package Generator
# Creates a distributable zip file with all necessary files

PACKAGE_VERSION="${1:-0.1.0}"
PACKAGE_NAME="STAR-RSeQC-${PACKAGE_VERSION}"
TEMP_DIR="/tmp/${PACKAGE_NAME}"
OUTPUT_DIR="$(pwd)"

echo "=== Generating STAR-RSeQC Distribution Package ==="
echo "Version: $PACKAGE_VERSION"
echo "Output: $OUTPUT_DIR/${PACKAGE_NAME}.zip"
echo

# ─── Cleanup ─────────────────────────────────────────────────────────────────

if [[ -d "$TEMP_DIR" ]]; then
    rm -rf "$TEMP_DIR"
fi
mkdir -p "$TEMP_DIR"

# ─── Build Release Binary ────────────────────────────────────────────────────

echo "[1/7] Building release binary..."
cargo build --release --locked
cp target/release/star-rseqc "$TEMP_DIR/star-rseqc"
strip "$TEMP_DIR/star-rseqc" || true
echo "✓ Binary built and copied"

# ─── Copy Scripts ────────────────────────────────────────────────────────────

echo "[2/7] Copying installation scripts..."
cp setup.sh "$TEMP_DIR/setup.sh"
chmod +x "$TEMP_DIR/setup.sh"
echo "✓ Scripts copied"

# ─── Copy Docker Files ───────────────────────────────────────────────────────

echo "[3/7] Copying Docker files..."
mkdir -p "$TEMP_DIR/docker"
cp Dockerfile "$TEMP_DIR/docker/"
cp .dockerignore "$TEMP_DIR/docker/"
cp docker-compose.yml "$TEMP_DIR/docker/"
echo "✓ Docker files copied"

# ─── Copy Documentation ─────────────────────────────────────────────────────

echo "[4/7] Copying documentation..."
mkdir -p "$TEMP_DIR/docs"
cp README.md "$TEMP_DIR/docs/"
cat > "$TEMP_DIR/docs/QUICK_START.md" << 'EOF'
# Quick Start

## Three Installation Options

### Option 1: Docker (Fastest, Recommended)
```bash
cd STAR-RSeQC-*
./setup.sh
# Select: 1 (Docker)
```

### Option 2: Automatic Setup (Mamba/Conda)
```bash
cd STAR-RSeQC-*
./setup.sh
# Select: 2 (Automatic)
```

### Option 3: Manual Installation
```bash
cd STAR-RSeQC-*
./setup.sh
# Select: 3 (Manual)
```

## Using the Binary Directly

If you've already installed dependencies:
```bash
./star-rseqc /path/to/fastq -o results
```

## Configuration

Edit `~/.config/star-rseqc/config.json`:
```json
{
  "genome_dir": "/path/to/star/index",
  "gtf": "/path/to/annotation.gtf",
  "star_env": "/path/to/conda/envs/star",
  "rseqc_env": "/path/to/conda/envs/rseqc",
  "samtools": "/path/to/samtools"
}
```

## Docker Usage

Build the image:
```bash
cd docker
docker build -t star-rseqc:latest .
```

Run analysis:
```bash
docker run --rm \
  -v /path/to/fastq:/data/input:ro \
  -v /path/to/output:/data/output \
  star-rseqc:latest \
  star-rseqc /data/input -o /data/output
```
EOF

cat > "$TEMP_DIR/docs/TROUBLESHOOTING.md" << 'EOF'
# Troubleshooting

## Docker Not Installed
Run the setup script and choose "Auto-install Docker" when prompted.

## Mamba/Conda Not Found
Run the setup script and choose "Automatic" for mamba installation.

## Binary Not in PATH
Add to your `.bashrc` or `.zshrc`:
```bash
export PATH="$HOME/.local/bin:$PATH"
```

## Config File Not Found
Create `~/.config/star-rseqc/config.json`:
```bash
mkdir -p ~/.config/star-rseqc
# Edit with your paths
nano ~/.config/star-rseqc/config.json
```

## STAR Environment Not Found
Verify conda environments:
```bash
conda env list
```

If missing, re-run setup.sh and choose automatic installation.

## Permission Denied
Make setup script executable:
```bash
chmod +x setup.sh
```

Make binary executable:
```bash
chmod +x star-rseqc
```

## Docker Permission Issues
Add your user to docker group:
```bash
sudo usermod -aG docker $USER
newgrp docker
```
EOF

echo "✓ Documentation copied"

# ─── Create Installation Instructions ────────────────────────────────────────

echo "[5/7] Creating installation instructions..."
cat > "$TEMP_DIR/INSTALL.txt" << 'EOF'
╔════════════════════════════════════════════════════════════════════════════╗
║                    STAR-RSeQC Installation Package                         ║
╚════════════════════════════════════════════════════════════════════════════╝

Thank you for downloading STAR-RSeQC!

═══════════════════════════════════════════════════════════════════════════════

QUICK START (Choose one):

  Option 1 - Docker (Recommended)
  ────────────────────────────────
  This requires Docker to be installed or auto-installable.

    $ ./setup.sh
    # Select: 1 (Docker)
    # Follow prompts to auto-install Docker if needed


  Option 2 - Automatic Mamba Setup
  ─────────────────────────────────
  Auto-installs mamba, conda environments, and all dependencies.

    $ ./setup.sh
    # Select: 2 (Automatic)
    # Follows all prompts


  Option 3 - Manual Installation
  ───────────────────────────────
  Follow web links and manual steps to install each component.

    $ ./setup.sh
    # Select: 3 (Manual)
    # Follow the provided links and instructions

═══════════════════════════════════════════════════════════════════════════════

WHAT'S INCLUDED:

  setup.sh              - Interactive installation script
  star-rseqc           - Pre-compiled binary
  docker/               - Dockerfile, docker-compose.yml, .dockerignore
  docs/                 - Full documentation and quick start guide

═══════════════════════════════════════════════════════════════════════════════

SYSTEM REQUIREMENTS:

  Minimum:
    - 8 GB RAM
    - 4 CPU cores
    - 50 GB disk space (for reference files)

  Recommended:
    - 32+ GB RAM
    - 16+ CPU cores
    - 200+ GB disk space

═══════════════════════════════════════════════════════════════════════════════

NEXT STEPS:

  1. Extract this archive
  2. Run: bash setup.sh
  3. Choose your installation method
  4. Follow the interactive prompts

═══════════════════════════════════════════════════════════════════════════════

DOCUMENTATION:

  Quick Start:        docs/QUICK_START.md
  Troubleshooting:    docs/TROUBLESHOOTING.md
  Full README:        docs/README.md

═══════════════════════════════════════════════════════════════════════════════

SUPPORT:

  GitHub:   https://github.com/Sandbox-commission/STAR-RSeQC
  Issues:   https://github.com/Sandbox-commission/STAR-RSeQC/issues

═══════════════════════════════════════════════════════════════════════════════
EOF

chmod +x "$TEMP_DIR/INSTALL.txt"
echo "✓ Installation instructions created"

# ─── Create Version File ──────────────────────────────────────────────────────

echo "[6/7] Creating version info..."
cat > "$TEMP_DIR/VERSION.txt" << EOF
STAR-RSeQC Distribution Package
Version: $PACKAGE_VERSION
Build Date: $(date -u)
Architecture: x86_64-linux-gnu
Build Type: Release

Binary: star-rseqc (stripped, optimized)

Key Components:
  - STAR 2.7.11b+ (via mamba)
  - RSeQC 5.0+ (via mamba)
  - samtools 1.15+ (via mamba)
  - Rust 1.70+ (build requirement only)

Installation Methods:
  1. Docker (requires Docker or auto-install)
  2. Automatic Mamba (auto-installs all dependencies)
  3. Manual (step-by-step web-based)
EOF
echo "✓ Version info created"

# ─── Create README for package ────────────────────────────────────────────────

echo "[7/7] Creating package README..."
cat > "$TEMP_DIR/README_PACKAGE.txt" << 'EOF'
STAR-RSeQC Distribution Contents
==================================

Directory Structure:

  star-rseqc              Pre-compiled binary
  setup.sh                Interactive installer (main entry point)
  docker/
    ├── Dockerfile        Container definition
    ├── docker-compose.yml Docker Compose configuration
    └── .dockerignore     Docker build exclusions
  docs/
    ├── README.md         Full documentation
    ├── QUICK_START.md    Getting started guide
    └── TROUBLESHOOTING.md Problem solver
  INSTALL.txt             This installation guide
  VERSION.txt             Version and build information

Getting Started:

  1. Make sure setup.sh is executable:
     chmod +x setup.sh

  2. Run the installer:
     ./setup.sh

  3. Choose your installation method:
     a) Docker (recommended, fastest)
     b) Automatic mamba/conda installation
     c) Manual installation with web links

  4. Follow the interactive prompts

Notes:

  - Installation typically takes 10-30 minutes (depending on method)
  - Docker method is fastest if Docker is already installed
  - Automatic mamba installation requires ~2-3 GB download
  - Manual installation provides maximum control

System Requirements:

  - Linux x86_64 or macOS (with Docker)
  - 8+ GB RAM minimum, 32+ GB recommended
  - Internet connection for downloads
  - ~50 GB disk space minimum

Troubleshooting:

  See docs/TROUBLESHOOTING.md for common issues and solutions.

Support:

  - GitHub Issues: https://github.com/Sandbox-commission/STAR-RSeQC/issues
  - Documentation: docs/README.md

EOF
echo "✓ Package README created"

# ─── Create ZIP Archive ──────────────────────────────────────────────────────

echo
echo "Creating archive: ${PACKAGE_NAME}.zip"
cd /tmp
zip -r -q "${OUTPUT_DIR}/${PACKAGE_NAME}.zip" "$PACKAGE_NAME"
rm -rf "$TEMP_DIR"

# ─── Summary ─────────────────────────────────────────────────────────────────

PACKAGE_SIZE=$(du -h "${OUTPUT_DIR}/${PACKAGE_NAME}.zip" | cut -f1)

echo
echo "═══════════════════════════════════════════════════════════════════════════"
echo "✓ Package created successfully!"
echo "═══════════════════════════════════════════════════════════════════════════"
echo
echo "Package Details:"
echo "  Name:     ${PACKAGE_NAME}.zip"
echo "  Location: ${OUTPUT_DIR}"
echo "  Size:     $PACKAGE_SIZE"
echo
echo "To use this package:"
echo "  1. Extract: unzip ${PACKAGE_NAME}.zip"
echo "  2. Enter:   cd ${PACKAGE_NAME}"
echo "  3. Run:     bash setup.sh"
echo
echo "Distribution:"
echo "  - Share the .zip file with users"
echo "  - Users extract and run setup.sh"
echo "  - No pre-installation required (setup handles dependencies)"
echo
