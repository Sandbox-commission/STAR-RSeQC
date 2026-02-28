use chrono::Local;
use crossterm::{
    cursor, event, execute,
    style::{self, Attribute, Color, Stylize},
    terminal::{self, ClearType},
};
use glob::glob;
use std::collections::HashMap;
use std::env;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Write as IoWrite};
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode, Stdio};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

static CANCELLED: AtomicBool = AtomicBool::new(false);

const REFRESH_INTERVAL: Duration = Duration::from_millis(100);

// ─── Configuration System ────────────────────────────────────────────────────
// Configuration resolution order (highest to lowest priority):
// 1. CLI flags (--genome-dir, --gtf, --star-env, etc.)
// 2. ~/.config/star-rseqc/config.json (user-local configuration file)
// 3. Auto-detection (find_conda_env, find_samtools)
// 4. Empty PathBuf (caught by validate_environment)

#[derive(serde::Deserialize, Default)]
struct FileConfig {
    genome_dir: Option<String>,
    gtf: Option<String>,
    #[allow(dead_code)]
    fasta: Option<String>,
    star_env: Option<String>,
    rseqc_env: Option<String>,
    samtools: Option<String>,
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
            Ok(content) => serde_json::from_str(&content).unwrap_or_default(),
            Err(_) => FileConfig::default(),
        },
    }
}

fn find_conda_env(name: &str) -> Option<PathBuf> {
    let candidates = ["miniforge3", "mambaforge", "miniconda3"];

    for cand in &candidates {
        if let Ok(home) = env::var("HOME") {
            let env_path = PathBuf::from(home).join(cand).join("envs").join(name);
            if env_path.exists() {
                return Some(env_path);
            }
        }
    }

    // Also check /opt/conda
    let opt_path = PathBuf::from("/opt/conda/envs").join(name);
    if opt_path.exists() {
        return Some(opt_path);
    }

    None
}

fn find_samtools(star_env: &Path) -> Option<PathBuf> {
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

fn usage() {
    eprintln!("star-rseqc v0.1.0");
    eprintln!("STAR 2-pass alignment + RSeQC quality control for paired-end RNA-seq");
    eprintln!();
    eprintln!("USAGE:");
    eprintln!("    star-rseqc <FASTQ_DIR> [OPTIONS]");
    eprintln!();
    eprintln!("DESCRIPTION:");
    eprintln!("    Discovers paired-end FASTQ samples (*_1P.fastq.gz / *_2P.fastq.gz)");
    eprintln!("    in the given directory and runs the following pipeline per sample:");
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
    eprintln!("    <FASTQ_DIR>               Directory containing *_1P.fastq.gz files");
    eprintln!("                              (paired R2 files must be *_2P.fastq.gz)");
    eprintln!();
    eprintln!("OPTIONS:");
    eprintln!("    -o, --output <DIR>        Output directory for all results");
    eprintln!("                              [default: star-rseqc-results]");
    eprintln!("    -j, --jobs <N>            Samples processed in parallel");
    eprintln!("                              [default: auto-detected from available RAM]");
    eprintln!("    -t, --threads <N>         Threads per STAR alignment job");
    eprintln!("                              [default: auto-detected from CPU count]");
    eprintln!("    --genome-dir <DIR>        STAR genome index directory");
    eprintln!("                              [set via --flag, ~/.config/star-rseqc/config.json, or auto-detected]");
    eprintln!("    --gtf <FILE>              GTF annotation file");
    eprintln!("                              [set via --flag, ~/.config/star-rseqc/config.json, or auto-detected]");
    eprintln!("    --bed <FILE>              Pre-built BED12 file for RSeQC");
    eprintln!("                              (auto-generated from GTF if omitted)");
    eprintln!("    --samtools <PATH>         Path to samtools binary");
    eprintln!("                              [set via --flag, ~/.config/star-rseqc/config.json, or auto-detected]");
    eprintln!("    --star-env <DIR>          STAR conda environment prefix (or 'auto' to search)");
    eprintln!("                              [set via --flag, ~/.config/star-rseqc/config.json, or auto-detected]");
    eprintln!("    --rseqc-env <DIR>         RSeQC conda environment prefix (or 'auto' to search)");
    eprintln!("                              [set via --flag, ~/.config/star-rseqc/config.json, or auto-detected]");
    eprintln!("    --bam-sort-ram <BYTES>    RAM limit for BAM sorting");
    eprintln!("                              [default: auto-detected from available RAM]");
    eprintln!("    --skip-qc                 Skip RSeQC steps (alignment only)");
    eprintln!("    --skip-alignment          Skip STAR (run QC on existing BAMs)");
    eprintln!("    --dry-run                 List samples without running anything");
    eprintln!("    -h, --help                Print this help message");
    eprintln!();
    eprintln!("FASTQ NAMING CONVENTION:");
    eprintln!("    Files must follow the pattern:");
    eprintln!("        <SAMPLE>_1P.fastq.gz   (read 1 / forward)");
    eprintln!("        <SAMPLE>_2P.fastq.gz   (read 2 / reverse)");
    eprintln!();
    eprintln!("    Sample name is everything before _1P / _2P, for example:");
    eprintln!("        103N_GBC_1P.fastq.gz  ->  sample = 103N_GBC");
    eprintln!("        50T_CRC_1P.fastq.gz   ->  sample = 50T_CRC");
    eprintln!();
    eprintln!("REFERENCE FILES:");
    eprintln!("    STAR index : [configured via --genome-dir or config file]");
    eprintln!("    GTF        : [configured via --gtf or config file]");
    eprintln!("    FASTA      : [used for building genome index if needed]");
    eprintln!();
    eprintln!("TOOL ENVIRONMENTS:");
    eprintln!("    STAR       : [auto-detected or set via --star-env]");
    eprintln!("    samtools   : [auto-detected or set via --samtools]");
    eprintln!("    RSeQC      : [auto-detected or set via --rseqc-env]");
    eprintln!();
    eprintln!("OUTPUT STRUCTURE:");
    eprintln!("    <output>/");
    eprintln!("      star/                          STAR alignment output per sample");
    eprintln!("        <sample>_Aligned.sortedByCoord.out.bam");
    eprintln!("        <sample>_Aligned.toTranscriptome.out.bam");
    eprintln!("        <sample>_ReadsPerGene.out.tab");
    eprintln!("        <sample>_Chimeric.out.junction");
    eprintln!("        <sample>_Log.final.out");
    eprintln!("      qc/                            RSeQC quality control output");
    eprintln!("        <sample>.strand.txt");
    eprintln!("        <sample>.geneBodyCoverage.txt");
    eprintln!("        <sample>.read_distribution.txt");
    eprintln!("      logs/                          Per-sample STAR log files");
    eprintln!("        <sample>.star.log");
    eprintln!("      annotation.bed12               Auto-generated BED12 (cached)");
    eprintln!("      pipeline_summary.json           JSON summary of all results");
    eprintln!("      pipeline_summary.tsv            TSV summary of all results");
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
    eprintln!("EXAMPLES:");
    eprintln!("    # Run on current directory (fastq.gz files here)");
    eprintln!("    star-rseqc ./");
    eprintln!();
    eprintln!("    # Run on a specific FASTQ directory");
    eprintln!("    star-rseqc /path/to/Paired/");
    eprintln!();
    eprintln!("    # Custom output and parallelism");
    eprintln!("    star-rseqc ./  -o my-results  -j 4  -t 8");
    eprintln!();
    eprintln!("    # Alignment only (skip QC)");
    eprintln!("    star-rseqc ./  --skip-qc");
    eprintln!();
    eprintln!("    # QC only on existing BAMs");
    eprintln!("    star-rseqc ./  --skip-alignment  -o existing-results/");
    eprintln!();
    eprintln!("    # Dry run to check sample discovery");
    eprintln!("    star-rseqc ./  --dry-run");
    eprintln!();
    eprintln!("    # Resume after interruption (just re-run)");
    eprintln!("    star-rseqc ./");
    eprintln!();
    eprintln!("NOTE:");
    eprintln!("    Resources (-j, -t, --bam-sort-ram) are auto-detected from system RAM");
    eprintln!("    and CPU count at startup. Each STAR job uses ~32 GB RAM for the genome.");
    eprintln!("    Override with explicit flags: -j 2 -t 16 --bam-sort-ram 4000000000");
    eprintln!("    Press Ctrl+C to gracefully cancel (waits for running jobs).");
}

// ─── System resource detection ───────────────────────────────────────────────

fn read_available_ram() -> u64 {
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

fn detect_system_resources() -> (u64, usize) {
    let ram = read_available_ram();
    let cpus = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    (ram, cpus)
}

fn auto_config_resources(available_ram: u64, total_cpus: usize) -> (usize, usize, u64) {
    const PER_JOB_RAM: u64 = 6_000_000_000; // 6 GB per job max
    const OS_BUFFER: u64 = 2_000_000_000; // reserve for OS
    let usable = available_ram.saturating_sub(OS_BUFFER);
    // Each job uses 6 GB total
    let jobs = ((usable / PER_JOB_RAM) as usize).max(1);
    let threads = (total_cpus / jobs).max(1);
    // BAM sort RAM capped at 6 GB per job
    let bam_sort = PER_JOB_RAM.min(usable / jobs as u64);
    (jobs, threads, bam_sort)
}

// ─── Config & Args ───────────────────────────────────────────────────────────

struct Config {
    fastq_dir: PathBuf,
    output_dir: PathBuf,
    genome_dir: PathBuf,
    gtf: PathBuf,
    bed: Option<PathBuf>,
    star_env: PathBuf,
    rseqc_env: PathBuf,
    samtools: PathBuf,
    threads_per_sample: usize,
    parallel_jobs: usize,
    bam_sort_ram: u64,
    skip_qc: bool,
    skip_alignment: bool,
    dry_run: bool,
    resources_auto: bool,
}

fn parse_args() -> Option<Config> {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        usage();
        return None;
    }

    // Load file configuration (layer 2 in priority order)
    let file_cfg = load_file_config();

    let mut fastq_dir: Option<PathBuf> = None;
    let mut output_dir = PathBuf::from("star-rseqc-results");
    // Initialize with file config layer (layer 2), will be overridden by CLI flags (layer 1)
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

    let mut i = 1;

    // Helper: consume the next argument as a value for a flag
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
    // For samtools, we need a valid star_env to search within it
    if samtools.is_none() {
        if let Some(ref star) = star_env {
            samtools = find_samtools(star);
        }
    }

    // Resolve resource parameters — auto-detect anything not explicitly set
    let resources_auto =
        parallel_jobs.is_none() || threads_per_sample.is_none() || bam_sort_ram.is_none();
    let (parallel_jobs, threads_per_sample, bam_sort_ram) = if resources_auto {
        let (avail_ram, total_cpus) = detect_system_resources();
        let (aj, at, ab) = auto_config_resources(avail_ram, total_cpus);
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
    })
}

fn is_cancelled() -> bool {
    CANCELLED.load(Ordering::Relaxed)
}

// ─── Sample types ────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
struct Sample {
    name: String,
    r1: PathBuf,
    r2: PathBuf,
}

#[derive(Clone)]
struct JobSlot {
    sample: String,
    step: String,
    started: Instant,
}

// ─── Progress state (msi-calc style) ─────────────────────────────────────────

struct ProgressState {
    total: usize,
    completed: AtomicUsize,
    skipped: AtomicUsize,
    failed: AtomicUsize,
    active_jobs: Mutex<Vec<Option<JobSlot>>>,
    recent_events: Mutex<Vec<String>>,
    phase_label: Mutex<String>,
    start_time: Instant,
    completed_durations: Mutex<Vec<f64>>,
}

impl ProgressState {
    fn new(total: usize, parallel_jobs: usize, phase: &str) -> Self {
        let slots = (0..parallel_jobs).map(|_| None).collect();
        Self {
            total,
            completed: AtomicUsize::new(0),
            skipped: AtomicUsize::new(0),
            failed: AtomicUsize::new(0),
            active_jobs: Mutex::new(slots),
            recent_events: Mutex::new(Vec::new()),
            phase_label: Mutex::new(phase.to_string()),
            start_time: Instant::now(),
            completed_durations: Mutex::new(Vec::new()),
        }
    }

    fn set_active(&self, slot: usize, sample: &str, step: &str) {
        if let Ok(mut jobs) = self.active_jobs.lock() {
            if slot < jobs.len() {
                jobs[slot] = Some(JobSlot {
                    sample: sample.to_string(),
                    step: step.to_string(),
                    started: Instant::now(),
                });
            }
        }
    }

    fn update_step(&self, slot: usize, step: &str) {
        if let Ok(mut jobs) = self.active_jobs.lock() {
            if slot < jobs.len() {
                if let Some(ref mut job) = jobs[slot] {
                    job.step = step.to_string();
                }
            }
        }
    }

    fn clear_slot(&self, slot: usize) {
        if let Ok(mut jobs) = self.active_jobs.lock() {
            if slot < jobs.len() {
                jobs[slot] = None;
            }
        }
    }

    fn add_event(&self, msg: String) {
        if let Ok(mut events) = self.recent_events.lock() {
            events.push(msg);
            if events.len() > 100 {
                events.remove(0);
            }
        }
    }

    fn done_count(&self) -> usize {
        self.completed.load(Ordering::Relaxed)
            + self.skipped.load(Ordering::Relaxed)
            + self.failed.load(Ordering::Relaxed)
    }

    fn record_duration(&self, secs: f64) {
        if let Ok(mut durations) = self.completed_durations.lock() {
            durations.push(secs);
        }
    }

    fn avg_duration(&self) -> f64 {
        self.completed_durations
            .lock()
            .ok()
            .and_then(|d| {
                if d.is_empty() {
                    None
                } else {
                    Some(d.iter().sum::<f64>() / d.len() as f64)
                }
            })
            .unwrap_or(0.0)
    }

    fn phase(&self) -> String {
        self.phase_label
            .lock()
            .map(|p| p.clone())
            .unwrap_or_default()
    }
}

// ─── SHA256 checkpoint system ────────────────────────────────────────────────
//
// On completion, a SHA256 digest is computed over the key output files for each
// sample (STAR Log.final.out + QC text files). These are small (KB-sized) so
// hashing is instant. The digest is stored in .checkpoints/<sample>.sha256.
//
// On resume, the digest is recomputed from the output files on disk. If any
// output was deleted, truncated, or corrupted, the hash won't match and the
// sample is automatically re-processed.
//
// This avoids hashing multi-GB FASTQ inputs while still providing cryptographic
// integrity verification of the pipeline results.

fn checkpoint_dir(output_dir: &Path) -> PathBuf {
    output_dir.join(".checkpoints")
}

/// SHA256 hash a single file (streamed, 64 KB chunks).
fn sha256_file(path: &Path) -> Result<Vec<u8>, String> {
    use sha2::{Digest, Sha256};
    use std::io::Read;

    let mut hasher = Sha256::new();
    let mut file = File::open(path)
        .map_err(|e| format!("Cannot open {} for hashing: {}", path.display(), e))?;
    let mut buf = [0u8; 65536];
    loop {
        let n = file
            .read(&mut buf)
            .map_err(|e| format!("Read error hashing {}: {}", path.display(), e))?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    Ok(hasher.finalize().to_vec())
}

/// Hash a list of files into a single SHA256 digest.
/// Each file's name and contents (or a __MISSING__ sentinel) are fed
/// into one streaming hasher, so the digest changes if any file is
/// added, removed, or modified.
fn sha256_file_list(files: &[PathBuf]) -> String {
    use sha2::{Digest, Sha256};

    let mut hasher = Sha256::new();
    for path in files {
        hasher.update(
            path.file_name()
                .unwrap_or_default()
                .to_string_lossy()
                .as_bytes(),
        );
        if path.exists() {
            if let Ok(bytes) = sha256_file(path) {
                hasher.update(&bytes);
            }
        } else {
            hasher.update(b"__MISSING__");
        }
    }
    format!("{:x}", hasher.finalize())
}

/// Directory-wise SHA256 digests for a sample's outputs.
struct SampleDigests {
    star: String,
    rseqc: String,
}

/// Compute separate SHA256 digests for STAR and RSeQC output directories.
///
/// STAR outputs (10 files — per STAR v2.7 manual):
///   {s}_Log.out, {s}_Log.progress.out, {s}_Log.final.out,
///   {s}_Aligned.sortedByCoord.out.bam, {s}_Aligned.sortedByCoord.out.bam.bai,
///   {s}_Aligned.toTranscriptome.out.bam, {s}_ReadsPerGene.out.tab,
///   {s}_SJ.out.tab, {s}_Chimeric.out.junction, {s}_Chimeric.out.sam
///
/// RSeQC outputs (6 files):
///   {s}.strand.txt, {s}.geneBodyCoverage.txt, {s}.geneBodyCoverage.r,
///   {s}.geneBodyCoverage.curves.pdf, {s}.geneBodyCoverage.heatMap.pdf,
///   {s}.read_distribution.txt
fn sha256_outputs(output_dir: &Path, sample_name: &str) -> SampleDigests {
    let star_dir = output_dir.join("star");
    let qc_dir = output_dir.join("qc");

    let star_files: Vec<PathBuf> = vec![
        star_dir.join(format!("{sample_name}_Log.out")),
        star_dir.join(format!("{sample_name}_Log.progress.out")),
        star_dir.join(format!("{sample_name}_Log.final.out")),
        star_dir.join(format!("{sample_name}_Aligned.sortedByCoord.out.bam")),
        star_dir.join(format!("{sample_name}_Aligned.sortedByCoord.out.bam.bai")),
        star_dir.join(format!("{sample_name}_Aligned.toTranscriptome.out.bam")),
        star_dir.join(format!("{sample_name}_ReadsPerGene.out.tab")),
        star_dir.join(format!("{sample_name}_SJ.out.tab")),
        star_dir.join(format!("{sample_name}_Chimeric.out.junction")),
        star_dir.join(format!("{sample_name}_Chimeric.out.sam")),
    ];

    let rseqc_files: Vec<PathBuf> = vec![
        qc_dir.join(format!("{sample_name}.strand.txt")),
        qc_dir.join(format!("{sample_name}.geneBodyCoverage.txt")),
        qc_dir.join(format!("{sample_name}.geneBodyCoverage.r")),
        qc_dir.join(format!("{sample_name}.geneBodyCoverage.curves.pdf")),
        qc_dir.join(format!("{sample_name}.geneBodyCoverage.heatMap.pdf")),
        qc_dir.join(format!("{sample_name}.read_distribution.txt")),
    ];

    SampleDigests {
        star: sha256_file_list(&star_files),
        rseqc: sha256_file_list(&rseqc_files),
    }
}

fn write_checkpoint(output_dir: &Path, name: &str, digests: &SampleDigests) {
    let dir = checkpoint_dir(output_dir);
    let _ = fs::create_dir_all(&dir);
    let path = dir.join(format!("{name}.sha256"));
    let _ = fs::write(
        &path,
        format!("star:{}\nrseqc:{}\n", digests.star, digests.rseqc),
    );
}

fn remove_checkpoint(output_dir: &Path, name: &str) {
    let path = checkpoint_dir(output_dir).join(format!("{name}.sha256"));
    let _ = fs::remove_file(&path);
}

/// Parse a checkpoint file into (star_digest, rseqc_digest).
fn parse_checkpoint(content: &str) -> Option<(String, String)> {
    let mut star = None;
    let mut rseqc = None;
    for line in content.lines() {
        let line = line.trim();
        if let Some(hex) = line.strip_prefix("star:") {
            star = Some(hex.to_string());
        } else if let Some(hex) = line.strip_prefix("rseqc:") {
            rseqc = Some(hex.to_string());
        }
    }
    Some((star?, rseqc?))
}

/// Check if a sample's outputs are complete and intact.
/// Returns:
///   SameHash        — both STAR and RSeQC digests match → skip
///   StarChanged     — only STAR outputs differ
///   RseqcChanged    — only RSeQC outputs differ
///   BothChanged     — both directories differ
///   NotDone         — no checkpoint → must process
enum ResumeStatus {
    SameHash,
    StarChanged {
        old: String,
        new: String,
    },
    RseqcChanged {
        old: String,
        new: String,
    },
    BothChanged {
        old_star: String,
        new_star: String,
        old_rseqc: String,
        new_rseqc: String,
    },
    NotDone,
}

fn check_resume(output_dir: &Path, sample_name: &str) -> ResumeStatus {
    let ckpt = checkpoint_dir(output_dir).join(format!("{sample_name}.sha256"));
    let content = match fs::read_to_string(&ckpt) {
        Ok(s) => s,
        Err(_) => return ResumeStatus::NotDone,
    };

    let (old_star, old_rseqc) = match parse_checkpoint(&content) {
        Some(pair) => pair,
        None => return ResumeStatus::NotDone, // malformed checkpoint
    };

    let current = sha256_outputs(output_dir, sample_name);
    let star_ok = old_star == current.star;
    let rseqc_ok = old_rseqc == current.rseqc;

    match (star_ok, rseqc_ok) {
        (true, true) => ResumeStatus::SameHash,
        (false, true) => ResumeStatus::StarChanged {
            old: old_star,
            new: current.star,
        },
        (true, false) => ResumeStatus::RseqcChanged {
            old: old_rseqc,
            new: current.rseqc,
        },
        (false, false) => ResumeStatus::BothChanged {
            old_star,
            new_star: current.star,
            old_rseqc,
            new_rseqc: current.rseqc,
        },
    }
}

// ─── Sample discovery ────────────────────────────────────────────────────────

fn discover_samples(fastq_dir: &Path) -> Vec<Sample> {
    let pattern = fastq_dir
        .join("*_1P.fastq.gz")
        .to_string_lossy()
        .to_string();

    let mut samples = Vec::new();
    let mut seen = HashMap::new();

    let entries: Vec<_> = match glob(&pattern) {
        Ok(paths) => paths.filter_map(|e| e.ok()).collect(),
        Err(_) => return samples,
    };

    for r1 in entries {
        let r1_name = r1.file_name().unwrap().to_string_lossy().to_string();

        let sample_name = match r1_name.strip_suffix("_1P.fastq.gz") {
            Some(n) => n.to_string(),
            None => continue,
        };

        let r2_name = format!("{}_2P.fastq.gz", sample_name);
        let r2 = r1.parent().unwrap().join(&r2_name);

        if !r2.exists() {
            eprintln!(
                "Warning: skipping {} — R2 not found ({})",
                sample_name,
                r2.display()
            );
            continue;
        }

        if seen.contains_key(&sample_name) {
            eprintln!("Warning: duplicate sample name {} — skipping", sample_name);
            continue;
        }
        seen.insert(sample_name.clone(), true);

        samples.push(Sample {
            name: sample_name,
            r1,
            r2,
        });
    }

    samples.sort_by(|a, b| a.name.cmp(&b.name));
    samples
}

// ─── GTF → BED12 conversion ─────────────────────────────────────────────────

fn extract_attribute(attrs: &str, key: &str) -> Option<String> {
    let search = format!("{} \"", key);
    if let Some(pos) = attrs.find(&search) {
        let start = pos + search.len();
        if let Some(end) = attrs[start..].find('"') {
            return Some(attrs[start..start + end].to_string());
        }
    }
    None
}

fn gtf_to_bed12(gtf_path: &Path, bed_path: &Path) -> Result<usize, String> {
    let gtf_file = File::open(gtf_path)
        .map_err(|e| format!("Cannot open GTF {}: {}", gtf_path.display(), e))?;
    let reader = BufReader::new(gtf_file);

    // Collect exons per transcript: (chrom, strand, Vec<(start, end)>)
    let mut transcripts: HashMap<String, (String, String, Vec<(u64, u64)>)> = HashMap::new();

    for line in reader.lines() {
        let line = line.map_err(|e| format!("Read error: {}", e))?;
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 || fields[2] != "exon" {
            continue;
        }

        let chrom = fields[0];
        let start: u64 = fields[3].parse::<u64>().unwrap_or(0).saturating_sub(1);
        let end: u64 = fields[4].parse().unwrap_or(0);
        let strand = fields[6];
        let attrs = fields[8];

        let transcript_id = match extract_attribute(attrs, "transcript_id") {
            Some(id) if !id.is_empty() => id,
            _ => continue,
        };

        transcripts
            .entry(transcript_id)
            .or_insert_with(|| (chrom.to_string(), strand.to_string(), Vec::new()))
            .2
            .push((start, end));
    }

    let out_file = File::create(bed_path)
        .map_err(|e| format!("Cannot create BED {}: {}", bed_path.display(), e))?;
    let mut writer = BufWriter::new(out_file);

    let mut count = 0usize;
    for (tx_id, (chrom, strand, ref mut exons)) in &mut transcripts {
        if exons.is_empty() {
            continue;
        }
        exons.sort_by_key(|e| e.0);

        let tx_start = exons[0].0;
        let tx_end = exons.last().unwrap().1;
        let block_count = exons.len();
        let block_sizes: Vec<String> = exons.iter().map(|e| (e.1 - e.0).to_string()).collect();
        let block_starts: Vec<String> =
            exons.iter().map(|e| (e.0 - tx_start).to_string()).collect();

        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t0\t{}\t{}\t{}\t0\t{}\t{}\t{}",
            chrom,
            tx_start,
            tx_end,
            tx_id,
            strand,
            tx_start,
            tx_end,
            block_count,
            block_sizes.join(","),
            block_starts.join(","),
        )
        .map_err(|e| format!("Write error: {}", e))?;
        count += 1;
    }

    if count == 0 {
        return Err("GTF→BED12 produced zero transcripts".to_string());
    }
    Ok(count)
}

// ─── Run command with cancellation ───────────────────────────────────────────

fn run_cancellable(mut cmd: Command) -> Result<bool, String> {
    let mut child = cmd.spawn().map_err(|e| format!("Failed to launch: {e}"))?;

    loop {
        match child.try_wait() {
            Ok(Some(status)) => return Ok(status.success()),
            Ok(None) => {
                if is_cancelled() {
                    let _ = child.kill();
                    let _ = child.wait();
                    return Err("Cancelled".to_string());
                }
                std::thread::sleep(Duration::from_millis(200));
            }
            Err(e) => return Err(format!("Wait error: {e}")),
        }
    }
}

fn make_log_stdio(log_dir: &Path, name: &str) -> (Stdio, Stdio) {
    let log_path = log_dir.join(format!("{name}.log"));
    match File::create(&log_path) {
        Ok(f) => {
            let f2 = f.try_clone().unwrap();
            (Stdio::from(f), Stdio::from(f2))
        }
        Err(_) => (Stdio::null(), Stdio::null()),
    }
}

// ─── Pipeline steps ──────────────────────────────────────────────────────────

fn process_sample(
    sample: &Sample,
    config: &Config,
    bed_path: &Path,
    state: &ProgressState,
    slot: usize,
) -> Result<(), String> {
    if is_cancelled() {
        return Err("Cancelled".to_string());
    }

    let star_dir = config.output_dir.join("star");
    let qc_dir = config.output_dir.join("qc");
    let log_dir = config.output_dir.join("logs");

    let bam_path = star_dir.join(format!("{}_Aligned.sortedByCoord.out.bam", sample.name));
    let out_prefix = format!("{}/{}_", star_dir.display(), sample.name);
    let job_start = Instant::now();

    // ── Step 1: STAR Alignment ──
    if !config.skip_alignment {
        state.set_active(slot, &sample.name, "STAR alignment");

        if bam_path.exists() {
            // Resume: BAM exists, skip alignment
            state.add_event(format!(
                "  RESUME  {} — BAM exists, skipping STAR",
                sample.name
            ));
        } else {
            let star_bin = config.star_env.join("bin/STAR");
            let (stdout_cfg, stderr_cfg) =
                make_log_stdio(&log_dir, &format!("{}.star", sample.name));

            let mut cmd = Command::new(&star_bin);
            cmd.args([
                "--runThreadN",
                &config.threads_per_sample.to_string(),
                "--genomeDir",
                config.genome_dir.to_str().unwrap(),
                "--readFilesIn",
                sample.r1.to_str().unwrap(),
                sample.r2.to_str().unwrap(),
                "--readFilesCommand",
                "zcat",
                "--outFileNamePrefix",
                &out_prefix,
                "--outSAMtype",
                "BAM",
                "SortedByCoordinate",
                "--twopassMode",
                "Basic",
                "--quantMode",
                "TranscriptomeSAM",
                "GeneCounts",
                "--outSAMstrandField",
                "intronMotif",
                "--chimSegmentMin",
                "15",
                "--chimJunctionOverhangMin",
                "15",
                "--chimScoreMin",
                "10",
                "--chimScoreDropMax",
                "30",
                "--chimScoreSeparation",
                "10",
                "--chimOutType",
                "Junctions",
                "SeparateSAMold",
                "--alignSJDBoverhangMin",
                "1",
                "--alignSJoverhangMin",
                "8",
                "--outFilterMismatchNoverReadLmax",
                "0.04",
                "--alignIntronMin",
                "20",
                "--alignIntronMax",
                "1000000",
                "--alignMatesGapMax",
                "1000000",
                "--limitBAMsortRAM",
                &config.bam_sort_ram.to_string(),
                "--sjdbGTFfile",
                config.gtf.to_str().unwrap(),
            ])
            .stdout(stdout_cfg)
            .stderr(stderr_cfg);

            match run_cancellable(cmd) {
                Ok(true) => {
                    state.add_event(format!("  DONE  {} — STAR alignment", sample.name));
                }
                Ok(false) => {
                    cleanup_partial_star(&star_dir, &sample.name);
                    state.add_event(format!(
                        "  FAIL  {} — STAR alignment error (cleaned)",
                        sample.name
                    ));
                    return Err(format!("{}: STAR failed", sample.name));
                }
                Err(e) => {
                    cleanup_partial_star(&star_dir, &sample.name);
                    return Err(e);
                }
            }
        }
    }

    if is_cancelled() {
        return Err("Cancelled".to_string());
    }

    // ── Verify BAM ──
    if !bam_path.exists() {
        state.add_event(format!(
            "  FAIL  {} — BAM not found after alignment",
            sample.name
        ));
        return Err(format!(
            "{}: BAM not found: {}",
            sample.name,
            bam_path.display()
        ));
    }

    // ── Step 2: samtools index ──
    state.update_step(slot, "samtools index");
    let bai_path = PathBuf::from(format!("{}.bai", bam_path.display()));
    if !bai_path.exists() {
        let mut cmd = Command::new(&config.samtools);
        cmd.args([
            "index",
            "-@",
            &config.threads_per_sample.to_string(),
            bam_path.to_str().unwrap(),
        ]);
        cmd.stdout(Stdio::null()).stderr(Stdio::null());

        match run_cancellable(cmd) {
            Ok(true) => {
                state.add_event(format!("  DONE  {} — samtools index", sample.name));
            }
            Ok(false) => {
                state.add_event(format!("  FAIL  {} — samtools index error", sample.name));
                return Err(format!("{}: samtools index failed", sample.name));
            }
            Err(e) => return Err(e),
        }
    }

    if is_cancelled() {
        return Err("Cancelled".to_string());
    }

    // ── Step 3: RSeQC ──
    if !config.skip_qc {
        let rseqc_python = config.rseqc_env.join("bin/python");

        // 3a: infer_experiment.py
        state.update_step(slot, "RSeQC: strandedness");
        let strand_out = qc_dir.join(format!("{}.strand.txt", sample.name));
        if !strand_out.exists() {
            let script = config.rseqc_env.join("bin/infer_experiment.py");
            let mut cmd = Command::new(&rseqc_python);
            cmd.args([
                script.to_str().unwrap(),
                "-i",
                bam_path.to_str().unwrap(),
                "-r",
                bed_path.to_str().unwrap(),
            ]);
            let output = cmd
                .output()
                .map_err(|e| format!("{}: infer_experiment.py spawn: {}", sample.name, e))?;
            if output.status.success() {
                fs::write(&strand_out, &output.stdout)
                    .map_err(|e| format!("Write {}: {}", strand_out.display(), e))?;
                state.add_event(format!("  DONE  {} — infer_experiment", sample.name));
            } else {
                let stderr = String::from_utf8_lossy(&output.stderr);
                state.add_event(format!("  FAIL  {} — infer_experiment", sample.name));
                // Non-fatal: continue with other QC steps
                eprintln!("{}: infer_experiment.py failed: {}", sample.name, stderr);
            }
        }

        if is_cancelled() {
            return Err("Cancelled".to_string());
        }

        // 3b: geneBody_coverage.py
        state.update_step(slot, "RSeQC: gene body coverage");
        let genebody_marker = qc_dir.join(format!("{}.geneBodyCoverage.txt", sample.name));
        if !genebody_marker.exists() {
            let script = config.rseqc_env.join("bin/geneBody_coverage.py");
            let genebody_prefix = format!("{}/{}", qc_dir.display(), sample.name);
            let mut cmd = Command::new(&rseqc_python);
            cmd.args([
                script.to_str().unwrap(),
                "-r",
                bed_path.to_str().unwrap(),
                "-i",
                bam_path.to_str().unwrap(),
                "-o",
                &genebody_prefix,
            ]);
            cmd.stdout(Stdio::null()).stderr(Stdio::null());

            match run_cancellable(cmd) {
                Ok(true) => {
                    state.add_event(format!("  DONE  {} — geneBody_coverage", sample.name));
                }
                Ok(false) => {
                    state.add_event(format!("  FAIL  {} — geneBody_coverage", sample.name));
                }
                Err(e) if e == "Cancelled" => return Err(e),
                Err(_) => {
                    state.add_event(format!("  FAIL  {} — geneBody_coverage", sample.name));
                }
            }
        }

        if is_cancelled() {
            return Err("Cancelled".to_string());
        }

        // 3c: read_distribution.py
        state.update_step(slot, "RSeQC: read distribution");
        let readdist_out = qc_dir.join(format!("{}.read_distribution.txt", sample.name));
        if !readdist_out.exists() {
            let script = config.rseqc_env.join("bin/read_distribution.py");
            let mut cmd = Command::new(&rseqc_python);
            cmd.args([
                script.to_str().unwrap(),
                "-i",
                bam_path.to_str().unwrap(),
                "-r",
                bed_path.to_str().unwrap(),
            ]);
            let output = cmd
                .output()
                .map_err(|e| format!("{}: read_distribution.py spawn: {}", sample.name, e))?;
            if output.status.success() {
                fs::write(&readdist_out, &output.stdout)
                    .map_err(|e| format!("Write {}: {}", readdist_out.display(), e))?;
                state.add_event(format!("  DONE  {} — read_distribution", sample.name));
            } else {
                state.add_event(format!("  FAIL  {} — read_distribution", sample.name));
            }
        }
    }

    // ── Done ──
    let dur = job_start.elapsed().as_secs_f64();
    state.record_duration(dur);
    Ok(())
}

// ─── Partial cleanup on failure ──────────────────────────────────────────────

#[allow(dead_code)]
fn cleanup_partial_star(star_dir: &Path, sample_name: &str) {
    let suffixes = [
        "_Aligned.sortedByCoord.out.bam",
        "_Aligned.sortedByCoord.out.bam.bai",
        "_Aligned.toTranscriptome.out.bam",
        "_ReadsPerGene.out.tab",
        "_Log.final.out",
        "_Log.out",
        "_Log.progress.out",
        "_SJ.out.tab",
        "_Chimeric.out.junction",
        "_Chimeric.out.sam",
    ];
    for suffix in &suffixes {
        let path = star_dir.join(format!("{sample_name}{suffix}"));
        if path.exists() {
            let _ = fs::remove_file(&path);
        }
    }
    // Remove the _STARgenome and _STARpass1 temp dirs
    for subdir in ["_STARgenome", "_STARpass1", "_STARtmp"] {
        let dir = star_dir.join(format!("{sample_name}{subdir}"));
        if dir.exists() {
            let _ = fs::remove_dir_all(&dir);
        }
    }
}

// ─── TUI rendering ──────────────────────────────────────────────────────────

fn fmt_duration(d: Duration) -> String {
    let secs = d.as_secs();
    let h = secs / 3600;
    let m = (secs % 3600) / 60;
    let s = secs % 60;
    if h > 0 {
        format!("{h:02}:{m:02}:{s:02}")
    } else {
        format!("{m:02}:{s:02}")
    }
}

fn fmt_secs(s: f64) -> String {
    fmt_duration(Duration::from_secs_f64(s))
}

fn render_screen(
    stdout: &mut io::Stdout,
    state: &ProgressState,
    parallel_jobs: usize,
    resumed: usize,
) {
    let (term_w, term_h) = terminal::size().unwrap_or((80, 24));
    let w = term_w as usize;
    let h = term_h as usize;

    let elapsed = state.start_time.elapsed();
    let done = state.done_count();
    let total = state.total;
    let completed = state.completed.load(Ordering::Relaxed);
    let skipped = state.skipped.load(Ordering::Relaxed);
    let failed = state.failed.load(Ordering::Relaxed);
    let remaining = total.saturating_sub(done);
    let phase = state.phase();
    let avg_dur = state.avg_duration();

    let pct = if total > 0 { done * 100 / total } else { 0 };
    let processed = completed;
    let speed = if elapsed.as_secs() > 0 && processed > 0 {
        processed as f64 / (elapsed.as_secs_f64() / 60.0)
    } else {
        0.0
    };
    let eta = if processed > 0 && remaining > 0 {
        Duration::from_secs_f64(avg_dur * remaining as f64)
    } else if done > 0 && remaining > 0 {
        let per = elapsed.as_secs_f64() / done as f64;
        Duration::from_secs_f64(per * remaining as f64)
    } else {
        Duration::ZERO
    };

    let bar_width = w.saturating_sub(2).min(60);
    let filled = if total > 0 {
        bar_width * done / total
    } else {
        0
    };
    let empty = bar_width.saturating_sub(filled);
    let bar_filled: String = "\u{2588}".repeat(filled);
    let bar_empty: String = "\u{2591}".repeat(empty);

    let active_snapshot: Vec<Option<JobSlot>> = state
        .active_jobs
        .lock()
        .map(|j| j.clone())
        .unwrap_or_default();
    let active_count = active_snapshot.iter().filter(|s| s.is_some()).count();

    let events: Vec<String> = state
        .recent_events
        .lock()
        .map(|e| e.clone())
        .unwrap_or_default();

    let _ = execute!(
        stdout,
        cursor::MoveTo(0, 0),
        terminal::Clear(ClearType::All)
    );

    let border_top = format!("\u{2554}{}\u{2557}", "\u{2550}".repeat(w.saturating_sub(2)));
    let border_mid = format!("\u{2560}{}\u{2563}", "\u{2550}".repeat(w.saturating_sub(2)));
    let border_bot = format!("\u{255A}{}\u{255D}", "\u{2550}".repeat(w.saturating_sub(2)));

    let mut row: u16 = 0;

    macro_rules! bline {
        ($color:expr, $text:expr) => {{
            let _ = execute!(
                stdout,
                cursor::MoveTo(0, row),
                style::SetForegroundColor(Color::Cyan)
            );
            let txt: String = $text;
            let txt_len = txt.len().min(w.saturating_sub(2));
            print!("\u{2551}");
            let _ = execute!(stdout, style::SetForegroundColor($color));
            print!("{}", &txt[..txt_len]);
            let _ = execute!(stdout, style::SetForegroundColor(Color::Cyan));
            print!(
                "{}\u{2551}",
                " ".repeat(w.saturating_sub(2).saturating_sub(txt_len))
            );
            row += 1;
        }};
    }

    macro_rules! separator {
        () => {{
            let _ = execute!(
                stdout,
                cursor::MoveTo(0, row),
                style::SetForegroundColor(Color::Cyan)
            );
            print!("{border_mid}");
            row += 1;
        }};
    }

    // Top border
    let _ = execute!(
        stdout,
        cursor::MoveTo(0, row),
        style::SetForegroundColor(Color::Cyan)
    );
    print!("{border_top}");
    row += 1;

    // Header
    let title = "STAR-RSeQC v0.1.0";
    let subtitle = "STAR 2-Pass Alignment + RSeQC Quality Control | Paired-End RNA-seq";
    let pad_title = (w.saturating_sub(2).saturating_sub(title.len())) / 2;
    let pad_sub = (w.saturating_sub(2).saturating_sub(subtitle.len())) / 2;

    let _ = execute!(stdout, cursor::MoveTo(0, row));
    print!(
        "\u{2551}{}{}{}\u{2551}",
        " ".repeat(pad_title),
        title.with(Color::White).attribute(Attribute::Bold),
        " ".repeat(
            w.saturating_sub(2)
                .saturating_sub(pad_title)
                .saturating_sub(title.len())
        )
    );
    row += 1;

    let _ = execute!(
        stdout,
        cursor::MoveTo(0, row),
        style::SetForegroundColor(Color::Cyan)
    );
    let sub_len = subtitle.len().min(w.saturating_sub(2));
    let sub_display = &subtitle[..sub_len];
    print!(
        "\u{2551}{}{}{}\u{2551}",
        " ".repeat(pad_sub.min(w.saturating_sub(2).saturating_sub(sub_len))),
        sub_display.with(Color::DarkGrey),
        " ".repeat(
            w.saturating_sub(2)
                .saturating_sub(pad_sub.min(w.saturating_sub(2).saturating_sub(sub_len)))
                .saturating_sub(sub_len)
        )
    );
    row += 1;

    separator!();

    // Phase indicator
    let phase_line = format!("  {phase}");
    bline!(Color::Magenta, phase_line);

    if resumed > 0 {
        let resume_line =
            format!("  Resumed: {resumed} sample(s) already completed from previous run");
        bline!(Color::DarkYellow, resume_line);
    }

    separator!();

    // Overall progress bar
    bline!(Color::White, "  OVERALL PROGRESS".to_string());

    let _ = execute!(
        stdout,
        cursor::MoveTo(0, row),
        style::SetForegroundColor(Color::Cyan)
    );
    print!("\u{2551}  ");
    let _ = execute!(stdout, style::SetForegroundColor(Color::Green));
    print!("{bar_filled}");
    let _ = execute!(stdout, style::SetForegroundColor(Color::DarkGrey));
    print!("{bar_empty}");
    let _ = execute!(stdout, style::SetForegroundColor(Color::White));
    let pct_str = format!(" {:>3}%", pct);
    print!("{pct_str}");
    let used = 2 + bar_width + pct_str.len();
    let _ = execute!(stdout, style::SetForegroundColor(Color::Cyan));
    print!(
        "{}\u{2551}",
        " ".repeat(w.saturating_sub(2).saturating_sub(used))
    );
    row += 1;

    let stats_line = format!(
        "  {}/{} done   Elapsed: {}   ETA: {}   Speed: {:.1}/min",
        done,
        total,
        fmt_duration(elapsed),
        fmt_duration(eta),
        speed
    );
    let stats_len = stats_line.len();
    let _ = execute!(
        stdout,
        cursor::MoveTo(0, row),
        style::SetForegroundColor(Color::Cyan)
    );
    print!(
        "\u{2551}{}{}\u{2551}",
        stats_line.with(Color::White),
        " ".repeat(w.saturating_sub(2).saturating_sub(stats_len))
    );
    row += 1;

    separator!();

    // Active jobs with per-sample progress bars
    let active_label = format!("  ACTIVE JOBS ({}/{})", active_count, parallel_jobs);
    let active_label_len = active_label.len();
    let _ = execute!(
        stdout,
        cursor::MoveTo(0, row),
        style::SetForegroundColor(Color::Cyan)
    );
    print!(
        "\u{2551}{}{}\u{2551}",
        active_label.with(Color::Yellow),
        " ".repeat(w.saturating_sub(2).saturating_sub(active_label_len))
    );
    row += 1;

    let active_jobs: Vec<(usize, &JobSlot)> = active_snapshot
        .iter()
        .enumerate()
        .filter_map(|(i, s)| s.as_ref().map(|j| (i, j)))
        .collect();

    let max_active_rows = (h.saturating_sub(row as usize + 12)) / 2;
    let spinner_chars = ['|', '/', '-', '\\'];
    let spin_idx = (elapsed.as_millis() / 250) as usize;
    let sample_bar_w = w.saturating_sub(8).min(40);

    if active_jobs.is_empty() {
        bline!(Color::DarkGrey, "  No active jobs".to_string());
    }

    for (shown, (i, job)) in active_jobs.iter().enumerate() {
        if shown >= max_active_rows {
            let hidden = active_jobs.len().saturating_sub(shown);
            if hidden > 0 {
                let more = format!("  ... and {hidden} more active");
                let more_len = more.len();
                let _ = execute!(
                    stdout,
                    cursor::MoveTo(0, row),
                    style::SetForegroundColor(Color::Cyan)
                );
                print!(
                    "\u{2551}{}{}\u{2551}",
                    more.with(Color::DarkGrey),
                    " ".repeat(w.saturating_sub(2).saturating_sub(more_len))
                );
                row += 1;
            }
            break;
        }

        let spin = spinner_chars[(spin_idx + i) % 4];
        let job_elapsed_secs = job.started.elapsed().as_secs_f64();
        let job_elapsed_str = fmt_secs(job_elapsed_secs);
        let max_name = w.saturating_sub(30).min(30);
        let name = if job.sample.len() > max_name {
            format!("{}...", &job.sample[..max_name.saturating_sub(3)])
        } else {
            job.sample.clone()
        };

        // Row 1: spinner + name + step + elapsed
        let eta_part = if avg_dur > 0.0 {
            format!("{} / ~{}", job_elapsed_str, fmt_secs(avg_dur))
        } else {
            job_elapsed_str
        };
        let step_label = format!("[{}]", job.step);
        let line = format!("  {spin} {name:<max_name$} {step_label} {eta_part}");
        let display = if line.len() > w.saturating_sub(2) {
            line[..w.saturating_sub(2)].to_string()
        } else {
            line
        };
        let _ = execute!(
            stdout,
            cursor::MoveTo(0, row),
            style::SetForegroundColor(Color::Cyan)
        );
        print!("\u{2551}");
        let _ = execute!(stdout, style::SetForegroundColor(Color::White));
        print!("{display}");
        let _ = execute!(stdout, style::SetForegroundColor(Color::Cyan));
        print!(
            "{}\u{2551}",
            " ".repeat(w.saturating_sub(2).saturating_sub(display.len()))
        );
        row += 1;

        // Row 2: per-sample progress bar
        let _ = execute!(
            stdout,
            cursor::MoveTo(0, row),
            style::SetForegroundColor(Color::Cyan)
        );
        print!("\u{2551}    ");
        if avg_dur > 0.0 {
            let frac = (job_elapsed_secs / avg_dur).min(1.0);
            let s_filled = (sample_bar_w as f64 * frac) as usize;
            let s_empty = sample_bar_w.saturating_sub(s_filled);
            let _ = execute!(stdout, style::SetForegroundColor(Color::Yellow));
            print!("{}", "\u{2588}".repeat(s_filled));
            let _ = execute!(stdout, style::SetForegroundColor(Color::DarkGrey));
            print!("{}", "\u{2591}".repeat(s_empty));
            let _ = execute!(stdout, style::SetForegroundColor(Color::White));
            let s_pct = format!(" {:>3}%", (frac * 100.0) as usize);
            print!("{s_pct}");
            let bar_used = 4 + sample_bar_w + s_pct.len();
            let _ = execute!(stdout, style::SetForegroundColor(Color::Cyan));
            print!(
                "{}\u{2551}",
                " ".repeat(w.saturating_sub(2).saturating_sub(bar_used))
            );
        } else {
            // Indeterminate: pulse animation
            let pulse_pos = (spin_idx + i * 3) % (sample_bar_w + 4);
            for p in 0..sample_bar_w {
                if p >= pulse_pos.saturating_sub(2) && p <= pulse_pos {
                    let _ = execute!(stdout, style::SetForegroundColor(Color::Yellow));
                    print!("\u{2588}");
                } else {
                    let _ = execute!(stdout, style::SetForegroundColor(Color::DarkGrey));
                    print!("\u{2591}");
                }
            }
            let bar_used = 4 + sample_bar_w;
            let _ = execute!(stdout, style::SetForegroundColor(Color::Cyan));
            print!(
                "{}\u{2551}",
                " ".repeat(w.saturating_sub(2).saturating_sub(bar_used))
            );
        }
        row += 1;
    }

    separator!();

    // Counters
    let counters = format!(
        "  Completed: {}   Skipped: {}   Failed: {}   Remaining: {}",
        completed, skipped, failed, remaining
    );
    let _ = execute!(
        stdout,
        cursor::MoveTo(0, row),
        style::SetForegroundColor(Color::Cyan)
    );
    print!("\u{2551}");
    let _ = execute!(stdout, style::SetForegroundColor(Color::Green));
    print!("  Completed: {completed}");
    let _ = execute!(stdout, style::SetForegroundColor(Color::Yellow));
    print!("   Skipped: {skipped}");
    let _ = execute!(stdout, style::SetForegroundColor(Color::Red));
    print!("   Failed: {failed}");
    let _ = execute!(stdout, style::SetForegroundColor(Color::White));
    print!("   Remaining: {remaining}");
    let used_len = counters.len();
    let _ = execute!(stdout, style::SetForegroundColor(Color::Cyan));
    print!(
        "{}\u{2551}",
        " ".repeat(w.saturating_sub(2).saturating_sub(used_len))
    );
    row += 1;

    separator!();

    // Recent activity
    let log_label = "  RECENT ACTIVITY";
    let _ = execute!(
        stdout,
        cursor::MoveTo(0, row),
        style::SetForegroundColor(Color::Cyan)
    );
    print!(
        "\u{2551}{}{}\u{2551}",
        log_label.with(Color::Magenta),
        " ".repeat(w.saturating_sub(2).saturating_sub(log_label.len()))
    );
    row += 1;

    let max_event_rows = h.saturating_sub(row as usize + 3);
    let start = events.len().saturating_sub(max_event_rows);
    for event_line in &events[start..] {
        let _ = execute!(
            stdout,
            cursor::MoveTo(0, row),
            style::SetForegroundColor(Color::Cyan)
        );
        let ev = if event_line.len() > w.saturating_sub(2) {
            event_line[..w.saturating_sub(2)].to_string()
        } else {
            event_line.clone()
        };
        print!("\u{2551}");
        if ev.contains("DONE") {
            let _ = execute!(stdout, style::SetForegroundColor(Color::Green));
        } else if ev.contains("SKIP") || ev.contains("RESUME") {
            let _ = execute!(stdout, style::SetForegroundColor(Color::Yellow));
        } else if ev.contains("FAIL") {
            let _ = execute!(stdout, style::SetForegroundColor(Color::Red));
        } else if ev.contains("STOP") {
            let _ = execute!(stdout, style::SetForegroundColor(Color::DarkRed));
        } else if ev.contains("INFO") {
            let _ = execute!(stdout, style::SetForegroundColor(Color::Cyan));
        } else {
            let _ = execute!(stdout, style::SetForegroundColor(Color::White));
        }
        print!("{ev}");
        let _ = execute!(stdout, style::SetForegroundColor(Color::Cyan));
        print!(
            "{}\u{2551}",
            " ".repeat(w.saturating_sub(2).saturating_sub(ev.len()))
        );
        row += 1;
    }

    // Fill remaining
    while (row as usize) < h.saturating_sub(2) {
        let _ = execute!(
            stdout,
            cursor::MoveTo(0, row),
            style::SetForegroundColor(Color::Cyan)
        );
        print!("\u{2551}{}\u{2551}", " ".repeat(w.saturating_sub(2)));
        row += 1;
    }

    // Footer
    let _ = execute!(
        stdout,
        cursor::MoveTo(0, row),
        style::SetForegroundColor(Color::Cyan)
    );
    let cancel_hint = if is_cancelled() {
        "  CANCELLING..."
    } else {
        "  Ctrl+C to cancel"
    };
    let timestamp = format!("Updated: {} ", Local::now().format("%H:%M:%S"));
    let footer_pad = w
        .saturating_sub(2)
        .saturating_sub(cancel_hint.len())
        .saturating_sub(timestamp.len());
    if is_cancelled() {
        print!(
            "\u{2551}{}{}{}\u{2551}",
            cancel_hint.with(Color::Red),
            " ".repeat(footer_pad),
            timestamp.with(Color::DarkGrey)
        );
    } else {
        print!(
            "\u{2551}{}{}{}\u{2551}",
            cancel_hint.with(Color::DarkGrey),
            " ".repeat(footer_pad),
            timestamp.with(Color::DarkGrey)
        );
    }
    row += 1;

    let _ = execute!(
        stdout,
        cursor::MoveTo(0, row),
        style::SetForegroundColor(Color::Cyan)
    );
    print!("{border_bot}");

    let _ = execute!(stdout, style::ResetColor);
    let _ = stdout.flush();
}

// ─── Display thread ──────────────────────────────────────────────────────────

struct DisplayThread {
    flag: Arc<AtomicBool>,
    handle: Option<std::thread::JoinHandle<()>>,
}

impl DisplayThread {
    fn start(state: Arc<ProgressState>, parallel_jobs: usize, resumed: usize) -> Self {
        let flag = Arc::new(AtomicBool::new(false));
        let display_flag = Arc::clone(&flag);

        let handle = std::thread::spawn(move || {
            let mut out = io::stdout();
            loop {
                render_screen(&mut out, &state, parallel_jobs, resumed);

                if event::poll(REFRESH_INTERVAL).unwrap_or(false) {
                    if let Ok(event::Event::Key(key)) = event::read() {
                        if key.code == event::KeyCode::Char('c')
                            && key.modifiers.contains(event::KeyModifiers::CONTROL)
                        {
                            CANCELLED.store(true, Ordering::Relaxed);
                        }
                    }
                }

                if display_flag.load(Ordering::Relaxed) || is_cancelled() {
                    render_screen(&mut out, &state, parallel_jobs, resumed);
                    break;
                }
            }
        });

        Self {
            flag,
            handle: Some(handle),
        }
    }

    fn stop(&mut self) {
        self.flag.store(true, Ordering::Relaxed);
        if let Some(h) = self.handle.take() {
            let _ = h.join();
        }
    }
}

// ─── Work queue ──────────────────────────────────────────────────────────────

fn run_work_queue<T, F>(
    items: &[T],
    parallel_jobs: usize,
    state: &Arc<ProgressState>,
    worker: F,
) -> Vec<String>
where
    T: Sync,
    F: Fn(&T, usize) -> Result<(), String> + Sync,
{
    let next_idx = AtomicUsize::new(0);
    let errors: Mutex<Vec<String>> = Mutex::new(Vec::new());

    std::thread::scope(|s| {
        for slot in 0..parallel_jobs {
            let next = &next_idx;
            let errs = &errors;
            let w = &worker;
            let st = &state;
            s.spawn(move || loop {
                if is_cancelled() {
                    break;
                }
                let idx = next.fetch_add(1, Ordering::Relaxed);
                if idx >= items.len() {
                    break;
                }
                if let Err(e) = w(&items[idx], slot) {
                    if e != "Cancelled" {
                        errs.lock().unwrap().push(e);
                    }
                }
                st.clear_slot(slot);
                std::thread::yield_now();
            });
        }
    });

    errors.into_inner().unwrap()
}

// ─── Environment validation ──────────────────────────────────────────────────

fn validate_environment(config: &Config) -> Result<(), String> {
    let star_bin = config.star_env.join("bin/STAR");
    if !star_bin.exists() {
        return Err(format!(
            "STAR binary not found: {}\nIs the --star-env path correct?",
            star_bin.display()
        ));
    }

    if !config.samtools.exists() {
        return Err(format!(
            "samtools not found: {}\nSpecify with --samtools",
            config.samtools.display()
        ));
    }

    for script in [
        "infer_experiment.py",
        "geneBody_coverage.py",
        "read_distribution.py",
    ] {
        let path = config.rseqc_env.join("bin").join(script);
        if !path.exists() {
            return Err(format!(
                "RSeQC script not found: {}\nIs the --rseqc-env path correct?",
                path.display()
            ));
        }
    }

    if !config.genome_dir.exists() {
        return Err(format!(
            "STAR genome dir not found: {}",
            config.genome_dir.display()
        ));
    }
    let genome_file = config.genome_dir.join("Genome");
    if !genome_file.exists() {
        return Err(format!(
            "STAR genome index incomplete (no Genome file in {})",
            config.genome_dir.display()
        ));
    }

    if !config.gtf.exists() {
        return Err(format!("GTF not found: {}", config.gtf.display()));
    }

    if !config.fastq_dir.exists() {
        return Err(format!(
            "FASTQ directory not found: {}",
            config.fastq_dir.display()
        ));
    }

    Ok(())
}

// ─── Main ────────────────────────────────────────────────────────────────────

fn main() -> ExitCode {
    let config = match parse_args() {
        Some(c) => c,
        None => return ExitCode::SUCCESS,
    };

    eprintln!("star-rseqc v0.1.0 | Run with -h or --help for usage information");
    eprintln!();

    // ── Validate environment ──
    if let Err(e) = validate_environment(&config) {
        eprintln!("Environment check failed:\n  {e}");
        return ExitCode::FAILURE;
    }
    eprintln!("Environment OK.");
    eprintln!(
        "Resources: {} job(s) x {} thread(s)/job, {:.1} GB BAM sort RAM{}",
        config.parallel_jobs,
        config.threads_per_sample,
        config.bam_sort_ram as f64 / 1e9,
        if config.resources_auto {
            " [auto-detected]"
        } else {
            " [manual]"
        }
    );

    // ── Discover samples ──
    let all_samples = discover_samples(&config.fastq_dir);
    if all_samples.is_empty() {
        eprintln!(
            "No paired-end samples found in {}\n\
             Expected files matching *_1P.fastq.gz with corresponding *_2P.fastq.gz",
            config.fastq_dir.display()
        );
        return ExitCode::FAILURE;
    }
    eprintln!("Discovered {} paired-end samples.", all_samples.len());

    // ── Create output structure ──
    for subdir in ["star", "qc", "logs"] {
        if let Err(e) = fs::create_dir_all(config.output_dir.join(subdir)) {
            eprintln!(
                "Cannot create {}/{}: {}",
                config.output_dir.display(),
                subdir,
                e
            );
            return ExitCode::FAILURE;
        }
    }
    let _ = fs::create_dir_all(checkpoint_dir(&config.output_dir));

    // ── BED12 file ──
    let bed_path = if let Some(ref bed) = config.bed {
        if !bed.exists() {
            eprintln!("BED file not found: {}", bed.display());
            return ExitCode::FAILURE;
        }
        bed.clone()
    } else {
        let auto_bed = config.output_dir.join("annotation.bed12");
        if auto_bed.exists() {
            eprintln!("Reusing cached BED12: {}", auto_bed.display());
        } else {
            eprintln!("Converting GTF → BED12...");
            match gtf_to_bed12(&config.gtf, &auto_bed) {
                Ok(n) => eprintln!("BED12: {} transcripts written.", n),
                Err(e) => {
                    eprintln!("GTF→BED12 failed: {e}");
                    return ExitCode::FAILURE;
                }
            }
        }
        auto_bed
    };

    // ── Resume detection (SHA256 verification) ──
    eprintln!("Checking resume status (SHA256 verification)...");
    let mut already_done: usize = 0;
    let mut output_changed: usize = 0;
    let mut to_process: Vec<&Sample> = Vec::new();

    for s in &all_samples {
        match check_resume(&config.output_dir, &s.name) {
            ResumeStatus::SameHash => {
                already_done += 1;
            }
            ResumeStatus::StarChanged { old, new } => {
                eprintln!(
                    "  STAR changed: {} (was {}..., now {}...)",
                    s.name,
                    &old[..12.min(old.len())],
                    &new[..12.min(new.len())]
                );
                output_changed += 1;
                to_process.push(s);
            }
            ResumeStatus::RseqcChanged { old, new } => {
                eprintln!(
                    "  RSeQC changed: {} (was {}..., now {}...)",
                    s.name,
                    &old[..12.min(old.len())],
                    &new[..12.min(new.len())]
                );
                output_changed += 1;
                to_process.push(s);
            }
            ResumeStatus::BothChanged {
                old_star,
                new_star,
                old_rseqc,
                new_rseqc,
            } => {
                eprintln!(
                    "  STAR+RSeQC changed: {} (star: {}…→{}…, rseqc: {}…→{}…)",
                    s.name,
                    &old_star[..12.min(old_star.len())],
                    &new_star[..12.min(new_star.len())],
                    &old_rseqc[..12.min(old_rseqc.len())],
                    &new_rseqc[..12.min(new_rseqc.len())]
                );
                output_changed += 1;
                to_process.push(s);
            }
            ResumeStatus::NotDone => {
                to_process.push(s);
            }
        }
    }

    if already_done > 0 {
        eprintln!(
            "Resuming: {already_done}/{} samples verified (SHA256 OK), {} to process.",
            all_samples.len(),
            to_process.len()
        );
    }
    if output_changed > 0 {
        eprintln!("  {output_changed} sample(s) have corrupted/changed outputs — will re-process.");
    }

    // ── Dry run ──
    if config.dry_run {
        println!();
        println!("Dry run — {} samples discovered:\n", all_samples.len());
        println!("  {:<25} {:<50} {}", "SAMPLE", "R1", "STATUS");
        println!("  {}", "-".repeat(90));
        for s in &all_samples {
            let status = match check_resume(&config.output_dir, &s.name) {
                ResumeStatus::SameHash => "DONE (SHA256 OK)",
                ResumeStatus::StarChanged { .. } => "STAR CHANGED (will re-run)",
                ResumeStatus::RseqcChanged { .. } => "RSeQC CHANGED (will re-run)",
                ResumeStatus::BothChanged { .. } => "STAR+RSeQC CHANGED (will re-run)",
                ResumeStatus::NotDone => "PENDING",
            };
            println!("  {:<25} {:<50} {}", s.name, s.r1.display(), status);
        }
        println!();
        println!("Resource plan:");
        let auto_tag = if config.resources_auto { " (auto)" } else { "" };
        println!(
            "  {} parallel job(s){} x {} thread(s)/job{} = {} total threads",
            config.parallel_jobs,
            auto_tag,
            config.threads_per_sample,
            auto_tag,
            config.parallel_jobs * config.threads_per_sample
        );
        println!(
            "  BAM sort RAM: {:.1} GB per job{}",
            config.bam_sort_ram as f64 / 1e9,
            auto_tag
        );
        println!("  Output: {}", config.output_dir.display());
        println!("  BED12:  {}", bed_path.display());
        return ExitCode::SUCCESS;
    }

    if to_process.is_empty() {
        eprintln!("All samples already completed. Nothing to do.");
        return ExitCode::SUCCESS;
    }

    eprintln!();
    std::thread::sleep(Duration::from_secs(2));

    // ── Enter TUI ──
    let mut stdout = io::stdout();
    let _ = execute!(stdout, terminal::EnterAlternateScreen, cursor::Hide);
    let _ = terminal::enable_raw_mode();

    let phase_label = format!(
        "STAR 2-pass + RSeQC ({} samples, {}x{}t{})",
        to_process.len(),
        config.parallel_jobs,
        config.threads_per_sample,
        if config.resources_auto { " [auto]" } else { "" }
    );
    let state = Arc::new(ProgressState::new(
        all_samples.len(),
        config.parallel_jobs,
        &phase_label,
    ));
    state.skipped.store(already_done, Ordering::Relaxed);
    if already_done > 0 {
        state.add_event(format!(
            "  RESUME  {already_done} sample(s) already completed"
        ));
    }

    let mut display = DisplayThread::start(Arc::clone(&state), config.parallel_jobs, already_done);

    let config_ref = &config;
    let bed_ref = &bed_path;

    let errors = run_work_queue(&to_process, config.parallel_jobs, &state, |sample, slot| {
        let result = process_sample(sample, config_ref, bed_ref, &state, slot);
        match &result {
            Ok(()) => {
                // Compute directory-wise SHA256 (STAR + RSeQC separately)
                let digests = sha256_outputs(&config_ref.output_dir, &sample.name);
                write_checkpoint(&config_ref.output_dir, &sample.name, &digests);
                state.completed.fetch_add(1, Ordering::Relaxed);
                state.add_event(format!(
                    "  DONE  {} — star:{}… rseqc:{}…",
                    sample.name,
                    &digests.star[..12],
                    &digests.rseqc[..12]
                ));
            }
            Err(e) if e == "Cancelled" => {
                state.add_event(format!("  STOP  {} — cancelled", sample.name));
            }
            Err(e) => {
                remove_checkpoint(&config_ref.output_dir, &sample.name);
                state.failed.fetch_add(1, Ordering::Relaxed);
                state.add_event(format!("  FAIL  {} — {}", sample.name, e));
            }
        }
        result
    });

    display.stop();

    // ── Leave TUI ──
    let _ = terminal::disable_raw_mode();
    let _ = execute!(stdout, cursor::Show, terminal::LeaveAlternateScreen);

    let was_cancelled = is_cancelled();
    let total = all_samples.len();
    let completed = state.completed.load(Ordering::Relaxed);
    let failed_count = state.failed.load(Ordering::Relaxed);
    let elapsed_str = fmt_duration(state.start_time.elapsed());

    // ── Write summary files ──
    write_summary_files(&config.output_dir, &all_samples);

    // ── Final summary ──
    println!();
    println!("  \u{2554}{}\u{2557}", "\u{2550}".repeat(52));
    if was_cancelled {
        println!("  \u{2551}         STAR-RSeQC  -  Cancelled by user           \u{2551}");
    } else {
        println!("  \u{2551}           STAR-RSeQC  -  Run Complete              \u{2551}");
    }
    println!("  \u{2560}{}\u{2563}", "\u{2550}".repeat(52));
    println!("  \u{2551}  Total samples:      {:<29}\u{2551}", total);
    println!("  \u{2551}  Completed (new):    {:<29}\u{2551}", completed);
    if already_done > 0 {
        println!(
            "  \u{2551}  Resumed (SHA256 OK): {:<28}\u{2551}",
            already_done
        );
    }
    if output_changed > 0 {
        println!(
            "  \u{2551}  Re-processed (corrupt):{:<26}\u{2551}",
            output_changed
        );
    }
    println!(
        "  \u{2551}  Failed:             {:<29}\u{2551}",
        failed_count
    );
    println!(
        "  \u{2551}  Elapsed:            {:<29}\u{2551}",
        elapsed_str
    );
    println!(
        "  \u{2551}  Threads/sample:     {:<29}\u{2551}",
        if config.resources_auto {
            format!("{} (auto)", config.threads_per_sample)
        } else {
            config.threads_per_sample.to_string()
        }
    );
    println!(
        "  \u{2551}  Parallel jobs:      {:<29}\u{2551}",
        if config.resources_auto {
            format!("{} (auto)", config.parallel_jobs)
        } else {
            config.parallel_jobs.to_string()
        }
    );
    println!("  \u{2560}{}\u{2563}", "\u{2550}".repeat(52));
    println!(
        "  \u{2551}  Output : {:<41}\u{2551}",
        config.output_dir.display()
    );
    println!(
        "  \u{2551}  BAMs   : {:<41}\u{2551}",
        config.output_dir.join("star").display()
    );
    println!(
        "  \u{2551}  QC     : {:<41}\u{2551}",
        config.output_dir.join("qc").display()
    );
    println!(
        "  \u{2551}  Logs   : {:<41}\u{2551}",
        config.output_dir.join("logs").display()
    );
    println!("  \u{255A}{}\u{255D}", "\u{2550}".repeat(52));
    println!();

    if !errors.is_empty() {
        println!("  Errors:");
        for err in &errors {
            println!("    - {err}");
        }
        println!();
    }

    if was_cancelled {
        println!("  Run was cancelled. Re-run the same command to resume.");
        println!();
    } else if !errors.is_empty() {
        println!("  Some samples failed. Re-run to retry failed samples.");
        println!();
    }

    if errors.is_empty() && !was_cancelled {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}

// ─── Summary files ───────────────────────────────────────────────────────────

fn write_summary_files(output_dir: &Path, samples: &[Sample]) {
    let star_dir = output_dir.join("star");
    let qc_dir = output_dir.join("qc");

    #[derive(serde::Serialize)]
    struct SummaryRow {
        sample: String,
        sha256: String,
        // STAR outputs
        log_final: bool,
        log_out: bool,
        log_progress: bool,
        bam_sorted: bool,
        bam_index: bool,
        bam_transcriptome: bool,
        gene_counts: bool,
        splice_junctions: bool,
        chimeric_junction: bool,
        chimeric_sam: bool,
        // RSeQC outputs
        strand_qc: bool,
        genebody_txt: bool,
        genebody_r: bool,
        genebody_curves_pdf: bool,
        genebody_heatmap_pdf: bool,
        readdist_qc: bool,
    }

    let ckpt_dir = checkpoint_dir(output_dir);
    let rows: Vec<SummaryRow> = samples
        .iter()
        .map(|s| {
            let n = &s.name;
            let sha_file = ckpt_dir.join(format!("{n}.sha256"));
            let (star_sha, rseqc_sha) = fs::read_to_string(&sha_file)
                .ok()
                .and_then(|c| parse_checkpoint(&c))
                .unwrap_or_default();
            let sha256 = format!(
                "star:{}|rseqc:{}",
                &star_sha[..16.min(star_sha.len())],
                &rseqc_sha[..16.min(rseqc_sha.len())]
            );
            SummaryRow {
                sample: n.clone(),
                sha256,
                // STAR
                log_final: star_dir.join(format!("{n}_Log.final.out")).exists(),
                log_out: star_dir.join(format!("{n}_Log.out")).exists(),
                log_progress: star_dir.join(format!("{n}_Log.progress.out")).exists(),
                bam_sorted: star_dir
                    .join(format!("{n}_Aligned.sortedByCoord.out.bam"))
                    .exists(),
                bam_index: star_dir
                    .join(format!("{n}_Aligned.sortedByCoord.out.bam.bai"))
                    .exists(),
                bam_transcriptome: star_dir
                    .join(format!("{n}_Aligned.toTranscriptome.out.bam"))
                    .exists(),
                gene_counts: star_dir.join(format!("{n}_ReadsPerGene.out.tab")).exists(),
                splice_junctions: star_dir.join(format!("{n}_SJ.out.tab")).exists(),
                chimeric_junction: star_dir.join(format!("{n}_Chimeric.out.junction")).exists(),
                chimeric_sam: star_dir.join(format!("{n}_Chimeric.out.sam")).exists(),
                // RSeQC
                strand_qc: qc_dir.join(format!("{n}.strand.txt")).exists(),
                genebody_txt: qc_dir.join(format!("{n}.geneBodyCoverage.txt")).exists(),
                genebody_r: qc_dir.join(format!("{n}.geneBodyCoverage.r")).exists(),
                genebody_curves_pdf: qc_dir
                    .join(format!("{n}.geneBodyCoverage.curves.pdf"))
                    .exists(),
                genebody_heatmap_pdf: qc_dir
                    .join(format!("{n}.geneBodyCoverage.heatMap.pdf"))
                    .exists(),
                readdist_qc: qc_dir.join(format!("{n}.read_distribution.txt")).exists(),
            }
        })
        .collect();

    // JSON
    if let Ok(json) = serde_json::to_string_pretty(&rows) {
        let _ = fs::write(output_dir.join("pipeline_summary.json"), &json);
    }

    // TSV (compact: group STAR and QC status)
    let mut tsv = String::from(
        "sample\tsha256\tlog_final\tbam_sorted\tbam_index\tbam_transcriptome\tgene_counts\t\
         splice_junctions\tchimeric_junction\tchimeric_sam\tstrand_qc\tgenebody_txt\t\
         genebody_r\tgenebody_curves_pdf\treaddist_qc\n",
    );
    for r in &rows {
        let ok = |b: bool| if b { "OK" } else { "MISSING" };
        tsv.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            r.sample,
            r.sha256,
            ok(r.log_final),
            ok(r.bam_sorted),
            ok(r.bam_index),
            ok(r.bam_transcriptome),
            ok(r.gene_counts),
            ok(r.splice_junctions),
            ok(r.chimeric_junction),
            ok(r.chimeric_sam),
            ok(r.strand_qc),
            ok(r.genebody_txt),
            ok(r.genebody_r),
            ok(r.genebody_curves_pdf),
            ok(r.readdist_qc),
        ));
    }
    let _ = fs::write(output_dir.join("pipeline_summary.tsv"), &tsv);
}
