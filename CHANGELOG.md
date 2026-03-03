# Changelog

All notable changes to STAR-RSeQC will be documented in this file.

## [0.2.0] - 2026-03-03

### Added
- **Flexible FASTQ naming**: `--r1-suffix` / `--r2-suffix` flags (default: `_1P`/`_2P`).
  Supports Illumina `_R1`/`_R2`, simple `_1`/`_2`, or any custom suffix.
- **`--star-extra-args`**: Pass additional STAR parameters as a quoted string,
  appended after the standard ENCODE-compliant parameters.
- **`--version` / `-V` flag**: Print version and exit.
- **Structured logging**: Uses `log` + `env_logger`. Set `RUST_LOG=debug` for verbose output.
- **CHANGELOG.md**: Version history tracking.
- **`docker-compose.override.yml.example`**: Separate debug configuration for Docker.

### Changed
- **Modularized source code**: Split 2,390-line `main.rs` into 7 modules
  (`config.rs`, `checkpoint.rs`, `sample.rs`, `gtf.rs`, `pipeline.rs`, `tui.rs`).
- **Improved error messages**: All tool launch failures now include the program name,
  path to relevant log files, and actionable hints (e.g., installation commands).
- **`docker-compose.yml`**: Removed `/bin/bash` command override so `docker-compose up`
  actually runs the pipeline. Resource limits now configurable via `CPU_LIMIT`,
  `MEM_LIMIT` environment variables.
- **`setup.sh`**: Fixed function-before-call ordering bug, replaced recursive `main`
  call with loop-based menu, changed `set -euo` to `set -eo` for interactive safety.

### Removed
- **`STAR-RSeQC-0.1.0.zip`**: Binary artifacts no longer tracked in git.
  Use GitHub Releases for distribution instead.
- **`create-package.sh`**: Redundant with GitHub Actions release workflow.

## [0.1.0] - 2025-01-01

### Added
- Initial release.
- STAR 2-pass alignment with ENCODE-compliant parameters.
- RSeQC quality control (infer_experiment, geneBody_coverage, read_distribution).
- SHA256 checkpoint system for resume support.
- Full-screen TUI with progress bars, ETA, and activity log.
- Dynamic CPU & RAM auto-detection.
- Three-tier installation system (Docker, Mamba, Manual).
- GTF-to-BED12 pure-Rust converter.
- JSON and TSV pipeline summary output.
- GitHub Actions CI/CD (build, lint, format, release).
