# STAR-RSeQC Multi-Tier Installation — Quick Reference

## TL;DR: Getting Started

```bash
# Extract package
unzip STAR-RSeQC-0.2.0.zip
cd STAR-RSeQC-0.2.0

# Run installer
bash setup.sh

# Choose path:
# [1] Docker (fastest, 5-15 min)
# [2] Mamba (auto-install, 20-30 min)  
# [3] Manual (full control, 30-60 min)
```

---

## Three Installation Paths

### 🐳 Path 1: Docker (Fastest)
- **Time**: 5-15 minutes
- **Dependencies**: Docker (auto-installable)
- **Best for**: Beginners, reproducible environments
- **What happens**:
  1. Check if Docker installed
  2. If no: Auto-install for Ubuntu/Debian/Fedora/RHEL/CentOS
  3. If yes: Build Docker image with STAR, RSeQC, samtools
  4. Create `docker-run.sh` script
  5. Done! Run: `./docker-run.sh /path/to/fastq`

### 🐍 Path 2: Automatic Mamba
- **Time**: 20-30 minutes  
- **Dependencies**: Bash + Internet (auto-installs mamba)
- **Best for**: Bioinformaticians, conda users
- **What happens**:
  1. Check for existing conda/mamba
  2. If no: Auto-download Miniforge, install to `~/miniforge3/`
  3. Create `star` environment (STAR 2.7.11b + samtools)
  4. Create `rseqc` environment (RSeQC + Python)
  5. Build binary: `cargo build --release`
  6. Write config: `~/.config/star-rseqc/config.json`
  7. Done! Run: `star-rseqc /path/to/fastq`

### 📋 Path 3: Manual
- **Time**: 30-60 minutes
- **Dependencies**: Bash + web links
- **Best for**: Advanced users, HPC clusters
- **What happens**:
  1. Script provides links to:
     - Miniforge: https://github.com/conda-forge/miniforge
     - Docker: https://docs.docker.com/engine/install/
  2. Step-by-step instructions for each component
  3. You install at your own pace
  4. User builds binary and creates config manually

---

## File Structure

```
STAR-RSeQC-0.2.0/
├── setup.sh                 ← Run this first
├── star-rseqc              ← Pre-compiled binary
├── docker/
│   ├── Dockerfile          ← Build Docker image
│   └── docker-compose.yml  ← Run with docker-compose
├── docs/
│   ├── README.md           ← Full documentation
│   ├── QUICK_START.md      ← Get started guide
│   └── TROUBLESHOOTING.md  ← Problem solving
├── INSTALL.txt             ← Installation guide
└── VERSION.txt             ← Version info
```

---

## Command Reference

### Start Installation
```bash
bash setup.sh
```

### Docker Path (after setup)
```bash
./docker-run.sh /path/to/fastq -o results
# Or with docker-compose:
docker-compose up
```

### Mamba Path (after setup)
```bash
export PATH="$HOME/.local/bin:$PATH"
star-rseqc /path/to/fastq -o results
```

### View Configuration
```bash
cat ~/.config/star-rseqc/config.json
```

### Manually Create Config
```bash
mkdir -p ~/.config/star-rseqc
cat > ~/.config/star-rseqc/config.json << 'END'
{
  "genome_dir": "/path/to/star/index",
  "gtf": "/path/to/annotation.gtf",
  "star_env": "/path/to/envs/star",
  "rseqc_env": "/path/to/envs/rseqc",
  "samtools": "/path/to/samtools"
}
END
```

---

## Comparison

| Aspect | Docker | Mamba | Manual |
|--------|--------|-------|--------|
| Speed | ⚡⚡⚡ (5-15m) | ⚡⚡ (20-30m) | ⚡ (30-60m) |
| Dependencies | Docker | Bash | Bash |
| Works Offline | ❌ | ❌ | ✅ |
| HPC Support | ❌ | ✅ | ✅ |
| Customization | Limited | Moderate | Full |
| Requires sudo | ✅ | ❌ | Maybe |

---

## Decision Guide

**Choose Docker if:**
- You want the fastest setup
- You want everything pre-configured
- You have Docker (or want to auto-install it)

**Choose Mamba if:**
- You're familiar with conda/bioconda
- You want moderate customization
- You don't want to use Docker

**Choose Manual if:**
- You need maximum control
- You're on an HPC cluster
- You want to use existing tools
- You prefer step-by-step instructions

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| Docker not installed | setup.sh will auto-install (Ubuntu/Debian/Fedora/RHEL) |
| Conda not found | setup.sh will auto-install Miniforge |
| Binary not in PATH | Add to `.bashrc`: `export PATH="$HOME/.local/bin:$PATH"` |
| Permission denied | `chmod +x setup.sh` |
| Can't find reference files | Edit `~/.config/star-rseqc/config.json` |
| Docker permission error | `sudo usermod -aG docker $USER` then `newgrp docker` |

See `docs/TROUBLESHOOTING.md` for more details.

---

## Support

- **Documentation**: See `docs/` directory
- **Installation Guide**: `INSTALLATION_GUIDE.md`
- **Quick Start**: `docs/QUICK_START.md`
- **Troubleshooting**: `docs/TROUBLESHOOTING.md`
- **GitHub Issues**: https://github.com/Sandbox-commission/STAR-RSeQC/issues

---

## System Requirements

| Path | Minimum | Recommended |
|------|---------|-------------|
| Docker | 8GB RAM, 4 CPUs | 32GB RAM, 16 CPUs |
| Mamba | 8GB RAM, 4 CPUs | 32GB RAM, 16 CPUs |
| Manual | 8GB RAM, 4 CPUs | 32GB RAM, 16 CPUs |

**Disk Space**: 50GB+ for reference files

---

## Package Creation (for maintainers)

```bash
./create-package.sh 0.2.0
# Creates: STAR-RSeQC-0.2.0.zip (480 KB)

# Share with users:
# 1. Upload ZIP to GitHub Releases
# 2. Users download and extract
# 3. Users run: bash setup.sh
# 4. Setup handles all dependencies
```

---

## Next Steps

1. ✅ Extract this package
2. ✅ Run `bash setup.sh`
3. ✅ Choose installation path (1, 2, or 3)
4. ✅ Follow prompts
5. ✅ Run STAR-RSeQC!

```bash
star-rseqc /path/to/fastq/directory
```
