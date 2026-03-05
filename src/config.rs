use std::env;
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::Command;

use log::{debug, info, warn};

pub(crate) const VERSION: &str = env!("CARGO_PKG_VERSION");

// ─── Configuration System ────────────────────────────────────────────────────
// Configuration resolution order (highest to lowest priority):
// 1. CLI flags (--genome-dir, --gtf, --star-env, etc.)
// 2. ~/.config/star-rseqc/config.json (user-local configuration file)
// 3. Auto-detection (find_conda_env, find_samtools)
// 4. Empty PathBuf (caught by validate_environment)

#[derive(serde::Deserialize, Default)]
pub(crate) struct FileConfig {
    pub genome_dir: Option<String>,
    pub gtf: Option<String>,
    #[allow(dead_code)]
    pub fasta: Option<String>,
    pub star_env: Option<String>,
    pub rseqc_env: Option<String>,
    pub samtools: Option<String>,
    pub r1_suffix: Option<String>,
    pub r2_suffix: Option<String>,
}

fn config_file_path() -> Option<PathBuf> {
    env::var_os("HOME").map(|home| {
        let home = PathBuf::from(home);
        home.join(".config/star-rseqc/config.json")
    })
}

fn load_file_config() -> FileConfig {
    match config_file_path() {
        None => FileConfig::default(),
        Some(path) => match fs::read_to_string(&path) {
            Ok(content) => {
                debug!("Loaded config from {}", path.display());
                match serde_json::from_str(&content) {
                    Ok(cfg) => cfg,
                    Err(e) => {
                        warn!(
                            "Failed to parse config file {}: {} — using defaults",
                            path.display(),
                            e
                        );
                        FileConfig::default()
                    }
                }
            }
            Err(_) => FileConfig::default(),
        },
    }
}

pub(crate) fn find_conda_env(name: &str) -> Option<PathBuf> {
    let candidates = ["miniforge3", "mambaforge", "miniconda3"];

    for cand in &candidates {
        if let Ok(home) = env::var("HOME") {
            let env_path = PathBuf::from(home).join(cand).join("envs").join(name);
            if env_path.exists() {
                debug!("Found conda env '{}' at {}", name, env_path.display());
                return Some(env_path);
            }
        }
    }

    // Also check /opt/conda
    let opt_path = PathBuf::from("/opt/conda/envs").join(name);
    if opt_path.exists() {
        debug!("Found conda env '{}' at {}", name, opt_path.display());
        return Some(opt_path);
    }

    None
}

pub(crate) fn find_samtools(star_env: &Path) -> Option<PathBuf> {
    let in_star_env = star_env.join("bin/samtools");
    if in_star_env.exists() {
        return Some(in_star_env);
    }

    // Check system PATH
    if let Ok(output) = Command::new("which").arg("samtools").output() {
        if output.status.success() {
            let path_str = String::from_utf8_lossy(&output.stdout);
            let path = path_str.trim();
            if !path.is_empty() {
                return Some(PathBuf::from(path));
            }
        }
    }

    None
}

// ─── Help ────────────────────────────────────────────────────────────────────

pub(crate) fn usage() {
    let (avail_ram, total_cpus) = detect_system_resources();
    let (auto_jobs, auto_threads, auto_bam_sort) =
        auto_config_resources(avail_ram, total_cpus, None);
    let ram_gb = avail_ram as f64 / 1e9;
    eprintln!("star-rseqc v{VERSION}");
    eprintln!("STAR 2-pass alignment + RSeQC quality control for paired-end RNA-seq");
    eprintln!();
    eprintln!("USAGE:");
    eprintln!("    star-rseqc <FASTQ_DIR> [OPTIONS]");
    eprintln!();
    eprintln!("DESCRIPTION:");
    eprintln!("    Discovers paired-end FASTQ samples in the given directory and runs");
    eprintln!("    the following pipeline per sample:");
    eprintln!();
    eprintln!("    1. STAR 2-pass alignment with chimeric junction detection");
    eprintln!("       - Produces coordinate-sorted BAM, transcriptome BAM, gene counts");
    eprintln!("       - Chimeric alignment for fusion detection");
    eprintln!("       - Passes --sjdbGTFfile for annotation-guided alignment");
    eprintln!();
    eprintln!("    2. samtools index on the sorted BAM");
    eprintln!();
    eprintln!("    3. RSeQC quality control:");
    eprintln!("       - infer_experiment.py  (library strandedness)");
    eprintln!("       - geneBody_coverage.py (5'-to-3' coverage uniformity)");
    eprintln!("       - read_distribution.py (genomic feature distribution)");
    eprintln!();
    eprintln!("    Samples are processed in parallel with a full-screen progress TUI.");
    eprintln!("    Resume-aware: re-run the same command to skip completed samples.");
    eprintln!();
    eprintln!("    If --genome-dir contains a FASTA file (.fa/.fna/.fasta) but no");
    eprintln!("    STAR index, the index is automatically generated before alignment.");
    eprintln!();
    eprintln!("    On first run, automatically converts the GTF annotation to BED12");
    eprintln!("    format required by RSeQC (cached for subsequent runs).");
    eprintln!();
    eprintln!("ARGUMENTS:");
    eprintln!("    <FASTQ_DIR>               Directory containing paired FASTQ files");
    eprintln!();
    eprintln!("OPTIONS:");
    eprintln!("    -o, --output <DIR>        Output directory for all results");
    eprintln!("                              [default: star-rseqc-results]");
    eprintln!("    -j, --jobs <N>            Samples processed in parallel");
    eprintln!(
        "                              [default: {auto_jobs} (auto: {total_cpus} CPUs, {ram_gb:.0} GB RAM)]"
    );
    eprintln!("    -t, --threads <N>         Threads per STAR alignment job");
    eprintln!(
        "                              [default: {auto_threads} (auto: {total_cpus} CPUs / {auto_jobs} jobs)]"
    );
    eprintln!("    --genome-dir <DIR>        STAR genome index directory");
    eprintln!("                              (auto-generates index if FASTA found but no index)");
    eprintln!("    --gtf <FILE>              GTF annotation file");
    eprintln!("    --bed <FILE>              Pre-built BED12 file for RSeQC");
    eprintln!("                              (auto-generated next to GTF file if omitted)");
    eprintln!("    --samtools <PATH>         Path to samtools binary");
    eprintln!("    --star-env <DIR>          STAR conda environment prefix (or 'auto')");
    eprintln!("    --rseqc-env <DIR>         RSeQC conda environment prefix (or 'auto')");
    eprintln!("    --bam-sort-ram <BYTES>    RAM limit for BAM sorting");
    eprintln!(
        "                              [default: {:.1} GB (auto)]",
        auto_bam_sort as f64 / 1e9
    );
    eprintln!("    --r1-suffix <SUFFIX>      Read 1 filename suffix before .fastq.gz");
    eprintln!("                              [default: _1P]  (e.g. _R1, _1)");
    eprintln!("    --r2-suffix <SUFFIX>      Read 2 filename suffix before .fastq.gz");
    eprintln!("                              [default: _2P]  (e.g. _R2, _2)");
    eprintln!("    --star-extra-args <ARGS>  Extra arguments passed to STAR (quoted string)");
    eprintln!("                              Appended after the standard ENCODE parameters");
    eprintln!("    --skip-qc                 Skip RSeQC steps (alignment only)");
    eprintln!("    --skip-alignment          Skip STAR (run QC on existing BAMs)");
    eprintln!("    --dry-run                 List samples without running anything");
    eprintln!("    -V, --version             Print version and exit");
    eprintln!("    -h, --help                Print this help message");
    eprintln!();
    eprintln!("FASTQ NAMING CONVENTION:");
    eprintln!("    By default, files must follow the pattern:");
    eprintln!("        <SAMPLE>_1P.fastq.gz   (read 1 / forward)");
    eprintln!("        <SAMPLE>_2P.fastq.gz   (read 2 / reverse)");
    eprintln!();
    eprintln!("    Use --r1-suffix / --r2-suffix for other conventions:");
    eprintln!("        --r1-suffix _R1 --r2-suffix _R2    (Illumina standard)");
    eprintln!("        --r1-suffix _1  --r2-suffix _2     (simple numbered)");
    eprintln!();
    eprintln!("STAR PARAMETERS:");
    eprintln!("    --twopassMode Basic               2-pass mapping for novel junctions");
    eprintln!("    --quantMode TranscriptomeSAM GeneCounts");
    eprintln!("    --outSAMstrandField intronMotif   Strand info for unstranded data");
    eprintln!("    --chimSegmentMin 15               Chimeric alignment for fusions");
    eprintln!("    --outFilterMismatchNoverReadLmax 0.04");
    eprintln!("    --alignIntronMax 1000000          Max intron length");
    eprintln!("    --alignMatesGapMax 1000000        Max mate pair gap");
    eprintln!();
    eprintln!("    Override or extend with --star-extra-args:");
    eprintln!(
        "        --star-extra-args \"--outFilterMultimapNmax 20 --winAnchorMultimapNmax 50\""
    );
    eprintln!();
    eprintln!("EXAMPLES:");
    eprintln!("    star-rseqc ./");
    eprintln!("    star-rseqc /path/to/Paired/");
    eprintln!("    star-rseqc ./  -o my-results  -j 4  -t 8");
    eprintln!("    star-rseqc ./  --skip-qc");
    eprintln!("    star-rseqc ./  --skip-alignment  -o existing-results/");
    eprintln!("    star-rseqc ./  --dry-run");
    eprintln!("    star-rseqc ./  --r1-suffix _R1 --r2-suffix _R2");
    eprintln!("    star-rseqc ./  --star-extra-args \"--outFilterMultimapNmax 20\"");
    eprintln!();
    eprintln!("NOTE:");
    eprintln!(
        "    Resources auto-detected: {total_cpus} CPUs, {ram_gb:.1} GB RAM available."
    );
    eprintln!(
        "    Auto plan: {auto_jobs} jobs x {auto_threads} threads = {} cores, {:.1} GB BAM sort/job.",
        auto_jobs * auto_threads,
        auto_bam_sort as f64 / 1e9
    );
    eprintln!("    STAR genome index RAM is detected from --genome-dir and excluded.");
    eprintln!("    Override with -j, -t, --bam-sort-ram if needed.");
    eprintln!("    Press Ctrl+C to gracefully cancel (waits for running jobs).");
}

// ─── System resource detection ───────────────────────────────────────────────

pub(crate) fn read_available_ram() -> u64 {
    let Ok(file) = std::fs::File::open("/proc/meminfo") else {
        return 0;
    };
    let reader = BufReader::new(file);
    for line in reader.lines().map_while(Result::ok) {
        if let Some(rest) = line.strip_prefix("MemAvailable:") {
            let kib: u64 = rest
                .split_whitespace()
                .next()
                .and_then(|s| s.parse().ok())
                .unwrap_or(0);
            return kib * 1024;
        }
    }
    0
}

pub(crate) fn detect_system_resources() -> (u64, usize) {
    let ram = read_available_ram();
    let cpus = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    (ram, cpus)
}

/// Estimate STAR genome index size by reading the `Genome` file.
/// Returns 0 if the file doesn't exist or can't be read.
fn estimate_genome_index_size(genome_dir: &Path) -> u64 {
    // The `Genome` file is the largest component; SA and SAindex add ~50% more.
    // Total shared memory footprint ≈ 2x the Genome file size.
    let genome_file = genome_dir.join("Genome");
    match fs::metadata(&genome_file) {
        Ok(m) => {
            let genome_bytes = m.len();
            debug!(
                "STAR Genome file: {:.1} GB -> estimated index footprint: {:.1} GB",
                genome_bytes as f64 / 1e9,
                (genome_bytes * 2) as f64 / 1e9
            );
            genome_bytes * 2 // Genome + SA + SAindex ≈ 2x Genome
        }
        Err(_) => 0,
    }
}

/// Compute optimal (jobs, threads, bam_sort_ram) from available system resources.
///
/// STAR resource model:
/// - **Shared memory**: STAR genome index is mmap'd and shared across all jobs.
///   For human GRCh38 this is ~27 GB. We detect this from the Genome file.
/// - **Per-job RAM**: Each STAR 2-pass job needs:
///   - BAM sort buffer (`--limitBAMsortRAM`, the main knob)
///   - Thread-local alignment buffers (~200 MB per thread)
///   - 2-pass intermediate data (~500 MB)
///   - RSeQC QC steps later (~500 MB)
/// - **OS reserve**: 2 GB for system, filesystem cache, etc.
///
/// Formula:
///   usable_ram = available_ram - genome_index - os_reserve
///   per_job_ram = bam_sort_ram + (threads * 200 MB) + 1 GB overhead
///   max_by_ram = usable_ram / per_job_ram
///   max_by_cpu = total_cpus / min_threads  (STAR needs ≥4 threads to be efficient)
///   jobs = min(max_by_ram, max_by_cpu)
///
/// Returns (parallel_jobs, threads_per_job, bam_sort_ram).
pub(crate) fn auto_config_resources(
    available_ram: u64,
    total_cpus: usize,
    genome_dir: Option<&Path>,
) -> (usize, usize, u64) {
    const OS_BUFFER: u64 = 2_000_000_000; // 2 GB for OS
    const DEFAULT_GENOME_INDEX: u64 = 32_000_000_000; // 32 GB fallback (human genome)
    const MIN_THREADS: usize = 4; // STAR needs ≥4 threads to be efficient
    const THREAD_BUFFER: u64 = 200_000_000; // ~200 MB per thread
    const JOB_OVERHEAD: u64 = 1_000_000_000; // 1 GB for 2-pass data + RSeQC
    const MIN_BAM_SORT: u64 = 2_000_000_000; // 2 GB minimum BAM sort
    const MAX_BAM_SORT: u64 = 10_000_000_000; // 10 GB maximum BAM sort

    // Detect genome index size, or use a conservative default
    let genome_index = genome_dir
        .map(estimate_genome_index_size)
        .filter(|&s| s > 0)
        .unwrap_or(DEFAULT_GENOME_INDEX);

    let usable_ram = available_ram
        .saturating_sub(genome_index)
        .saturating_sub(OS_BUFFER);

    if usable_ram == 0 {
        // Not enough RAM detected; return minimum viable config
        info!(
            "Low RAM detected ({:.1} GB available, {:.1} GB genome index). Using minimum config.",
            available_ram as f64 / 1e9,
            genome_index as f64 / 1e9
        );
        return (1, total_cpus.max(1), MIN_BAM_SORT);
    }

    // Start with max jobs by CPU (STAR needs MIN_THREADS per job)
    let max_by_cpu = (total_cpus / MIN_THREADS).max(1);

    // For each candidate job count, calculate per-job RAM budget
    // and find the highest job count that fits in RAM
    let mut best_jobs = 1;
    let mut best_threads = total_cpus;
    let mut best_bam_sort = MIN_BAM_SORT;

    for candidate_jobs in 1..=max_by_cpu {
        let threads = total_cpus / candidate_jobs;
        if threads < MIN_THREADS {
            break;
        }

        let per_job_overhead = (threads as u64) * THREAD_BUFFER + JOB_OVERHEAD;
        let ram_for_bam_sort = (usable_ram / candidate_jobs as u64).saturating_sub(per_job_overhead);

        if ram_for_bam_sort >= MIN_BAM_SORT {
            best_jobs = candidate_jobs;
            best_threads = threads;
            best_bam_sort = ram_for_bam_sort.min(MAX_BAM_SORT);
        }
    }

    debug!(
        "Resource plan: {:.0} GB available - {:.0} GB genome - {:.0} GB OS = {:.0} GB usable \
         -> {} jobs x {} threads, {:.1} GB BAM sort",
        available_ram as f64 / 1e9,
        genome_index as f64 / 1e9,
        OS_BUFFER as f64 / 1e9,
        usable_ram as f64 / 1e9,
        best_jobs,
        best_threads,
        best_bam_sort as f64 / 1e9,
    );

    (best_jobs, best_threads, best_bam_sort)
}

// ─── Config & Args ───────────────────────────────────────────────────────────

pub(crate) struct Config {
    pub fastq_dir: PathBuf,
    pub output_dir: PathBuf,
    pub genome_dir: PathBuf,
    pub gtf: PathBuf,
    pub bed: Option<PathBuf>,
    pub star_env: PathBuf,
    pub rseqc_env: PathBuf,
    pub samtools: PathBuf,
    pub threads_per_sample: usize,
    pub parallel_jobs: usize,
    pub bam_sort_ram: u64,
    pub skip_qc: bool,
    pub skip_alignment: bool,
    pub dry_run: bool,
    pub resources_auto: bool,
    pub r1_suffix: String,
    pub r2_suffix: String,
    pub star_extra_args: Vec<String>,
}

pub(crate) fn parse_args() -> Option<Config> {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        usage();
        return None;
    }

    // Load file configuration (layer 2 in priority order)
    let file_cfg = load_file_config();

    let mut fastq_dir: Option<PathBuf> = None;
    let mut output_dir = PathBuf::from("star-rseqc-results");
    let mut genome_dir: Option<PathBuf> = file_cfg.genome_dir.map(PathBuf::from);
    let mut gtf: Option<PathBuf> = file_cfg.gtf.map(PathBuf::from);
    let mut bed: Option<PathBuf> = None;
    let mut star_env: Option<PathBuf> = file_cfg.star_env.map(PathBuf::from);
    let mut rseqc_env: Option<PathBuf> = file_cfg.rseqc_env.map(PathBuf::from);
    let mut samtools: Option<PathBuf> = file_cfg.samtools.map(PathBuf::from);
    let mut threads_per_sample: Option<usize> = None;
    let mut parallel_jobs: Option<usize> = None;
    let mut bam_sort_ram: Option<u64> = None;
    let mut skip_qc = false;
    let mut skip_alignment = false;
    let mut dry_run = false;
    let mut r1_suffix = file_cfg.r1_suffix.unwrap_or_else(|| "_1P".to_string());
    let mut r2_suffix = file_cfg.r2_suffix.unwrap_or_else(|| "_2P".to_string());
    let mut star_extra_args: Vec<String> = Vec::new();

    let mut i = 1;

    macro_rules! next_val {
        ($flag:expr) => {{
            i += 1;
            if i >= args.len() {
                eprintln!("Error: {} requires a value.", $flag);
                return None;
            }
            &args[i]
        }};
    }
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "--help" => {
                usage();
                return None;
            }
            "-V" | "--version" => {
                eprintln!("star-rseqc v{VERSION}");
                return None;
            }
            "-o" | "--output" => {
                output_dir = PathBuf::from(next_val!("-o/--output"));
            }
            "-j" | "--jobs" => {
                let v = next_val!("-j/--jobs");
                parallel_jobs = Some(v.parse().unwrap_or_else(|_| {
                    eprintln!("Invalid value for --jobs: {v}");
                    std::process::exit(1);
                }));
            }
            "-t" | "--threads" => {
                let v = next_val!("-t/--threads");
                threads_per_sample = Some(v.parse().unwrap_or_else(|_| {
                    eprintln!("Invalid value for --threads: {v}");
                    std::process::exit(1);
                }));
            }
            "--genome-dir" => {
                genome_dir = Some(PathBuf::from(next_val!("--genome-dir")));
            }
            "--gtf" => {
                gtf = Some(PathBuf::from(next_val!("--gtf")));
            }
            "--bed" => {
                bed = Some(PathBuf::from(next_val!("--bed")));
            }
            "--star-env" => {
                let val = next_val!("--star-env");
                if val == "auto" {
                    star_env = find_conda_env("star");
                } else {
                    star_env = Some(PathBuf::from(val));
                }
            }
            "--rseqc-env" => {
                let val = next_val!("--rseqc-env");
                if val == "auto" {
                    rseqc_env = find_conda_env("rseqc");
                } else {
                    rseqc_env = Some(PathBuf::from(val));
                }
            }
            "--samtools" => {
                samtools = Some(PathBuf::from(next_val!("--samtools")));
            }
            "--bam-sort-ram" => {
                let v = next_val!("--bam-sort-ram");
                bam_sort_ram = Some(v.parse().unwrap_or_else(|_| {
                    eprintln!("Invalid value for --bam-sort-ram: {v}");
                    std::process::exit(1);
                }));
            }
            "--r1-suffix" => {
                r1_suffix = next_val!("--r1-suffix").to_string();
            }
            "--r2-suffix" => {
                r2_suffix = next_val!("--r2-suffix").to_string();
            }
            "--star-extra-args" => {
                let val = next_val!("--star-extra-args");
                star_extra_args = val.split_whitespace().map(|s| s.to_string()).collect();
            }
            "--skip-qc" => skip_qc = true,
            "--skip-alignment" => skip_alignment = true,
            "--dry-run" => dry_run = true,
            other => {
                if other.starts_with('-') {
                    eprintln!("Unknown option: {other}");
                    eprintln!("Run with -h for help.");
                    return None;
                }
                // Positional argument: FASTQ directory
                fastq_dir = Some(PathBuf::from(other));
            }
        }
        i += 1;
    }

    let fastq_dir = match fastq_dir {
        Some(d) => d,
        None => {
            eprintln!("Error: FASTQ_DIR argument is required.");
            eprintln!("Run with -h for help.");
            return None;
        }
    };

    // Apply layer 3 auto-detection: star_env, rseqc_env, samtools
    if star_env.is_none() {
        star_env = find_conda_env("star");
    }
    if rseqc_env.is_none() {
        rseqc_env = find_conda_env("rseqc");
    }
    if samtools.is_none() {
        if let Some(ref star) = star_env {
            samtools = find_samtools(star);
        }
    }

    // Resolve resource parameters
    let resources_auto =
        parallel_jobs.is_none() || threads_per_sample.is_none() || bam_sort_ram.is_none();
    let (parallel_jobs, threads_per_sample, bam_sort_ram) = if resources_auto {
        let (avail_ram, total_cpus) = detect_system_resources();
        let gd = genome_dir.as_deref().filter(|p| p.exists());
        let (aj, at, ab) = auto_config_resources(avail_ram, total_cpus, gd);
        info!(
            "Auto-detected resources: {:.1} GB RAM, {} CPUs -> {} jobs x {} threads, {:.1} GB BAM sort",
            avail_ram as f64 / 1e9,
            total_cpus,
            parallel_jobs.unwrap_or(aj),
            threads_per_sample.unwrap_or(at),
            bam_sort_ram.unwrap_or(ab) as f64 / 1e9,
        );
        (
            parallel_jobs.unwrap_or(aj),
            threads_per_sample.unwrap_or(at),
            bam_sort_ram.unwrap_or(ab),
        )
    } else {
        (
            parallel_jobs.unwrap(),
            threads_per_sample.unwrap(),
            bam_sort_ram.unwrap(),
        )
    };

    if parallel_jobs == 0 {
        eprintln!("Error: --jobs must be >= 1");
        return None;
    }
    if threads_per_sample == 0 {
        eprintln!("Error: --threads must be >= 1");
        return None;
    }

    if !star_extra_args.is_empty() {
        info!("Extra STAR args: {:?}", star_extra_args);
    }

    if r1_suffix.is_empty() {
        eprintln!("Error: --r1-suffix must be non-empty (e.g. _1P, _R1, _1)");
        return None;
    }
    if r2_suffix.is_empty() {
        eprintln!("Error: --r2-suffix must be non-empty (e.g. _2P, _R2, _2)");
        return None;
    }
    if r1_suffix.contains('/') || r1_suffix.contains('\\') {
        eprintln!("Error: --r1-suffix must not contain path separators");
        return None;
    }
    if r2_suffix.contains('/') || r2_suffix.contains('\\') {
        eprintln!("Error: --r2-suffix must not contain path separators");
        return None;
    }

    if r1_suffix != "_1P" || r2_suffix != "_2P" {
        info!(
            "Using custom FASTQ suffixes: R1='{}', R2='{}'",
            r1_suffix, r2_suffix
        );
    }

    Some(Config {
        fastq_dir,
        output_dir,
        genome_dir: genome_dir.unwrap_or_default(),
        gtf: gtf.unwrap_or_default(),
        bed,
        star_env: star_env.unwrap_or_default(),
        rseqc_env: rseqc_env.unwrap_or_default(),
        samtools: samtools.unwrap_or_default(),
        threads_per_sample,
        parallel_jobs,
        bam_sort_ram,
        skip_qc,
        skip_alignment,
        dry_run,
        resources_auto,
        r1_suffix,
        r2_suffix,
        star_extra_args,
    })
}
