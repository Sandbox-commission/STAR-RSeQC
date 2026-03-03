use std::env;
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::Command;

use log::{debug, info};

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
                serde_json::from_str(&content).unwrap_or_default()
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
    eprintln!("                              [default: auto-detected from available RAM]");
    eprintln!("    -t, --threads <N>         Threads per STAR alignment job");
    eprintln!("                              [default: auto-detected from CPU count]");
    eprintln!("    --genome-dir <DIR>        STAR genome index directory");
    eprintln!("    --gtf <FILE>              GTF annotation file");
    eprintln!("    --bed <FILE>              Pre-built BED12 file for RSeQC");
    eprintln!("                              (auto-generated from GTF if omitted)");
    eprintln!("    --samtools <PATH>         Path to samtools binary");
    eprintln!("    --star-env <DIR>          STAR conda environment prefix (or 'auto')");
    eprintln!("    --rseqc-env <DIR>         RSeQC conda environment prefix (or 'auto')");
    eprintln!("    --bam-sort-ram <BYTES>    RAM limit for BAM sorting");
    eprintln!("                              [default: auto-detected from available RAM]");
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
    eprintln!("        --star-extra-args \"--outFilterMultimapNmax 20 --winAnchorMultimapNmax 50\"");
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
    eprintln!("    Resources (-j, -t, --bam-sort-ram) are auto-detected from system RAM");
    eprintln!("    and CPU count at startup. Override with explicit flags.");
    eprintln!("    Press Ctrl+C to gracefully cancel (waits for running jobs).");
}

// ─── System resource detection ───────────────────────────────────────────────

pub(crate) fn read_available_ram() -> u64 {
    let Ok(file) = std::fs::File::open("/proc/meminfo") else {
        return 0;
    };
    let reader = BufReader::new(file);
    for line in reader.lines().flatten() {
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

pub(crate) fn auto_config_resources(available_ram: u64, total_cpus: usize) -> (usize, usize, u64) {
    const PER_JOB_RAM: u64 = 6_000_000_000; // 6 GB per job max
    const OS_BUFFER: u64 = 2_000_000_000; // reserve for OS
    let usable = available_ram.saturating_sub(OS_BUFFER);
    let jobs = ((usable / PER_JOB_RAM) as usize).max(1);
    let threads = (total_cpus / jobs).max(1);
    let bam_sort = PER_JOB_RAM.min(usable / jobs as u64);
    (jobs, threads, bam_sort)
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
        let (aj, at, ab) = auto_config_resources(avail_ram, total_cpus);
        info!(
            "Auto-detected resources: {} GB RAM, {} CPUs -> {} jobs x {} threads",
            avail_ram / 1_000_000_000,
            total_cpus,
            parallel_jobs.unwrap_or(aj),
            threads_per_sample.unwrap_or(at)
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
