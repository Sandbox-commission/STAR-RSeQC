use chrono::Local;
use crossterm::{
    cursor, event, execute,
    style::{self, Attribute, Color},
    terminal::{self, ClearType},
    tty::IsTty,
};
use glob::glob;
use std::collections::{HashMap, VecDeque};
use std::env;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write as IoWrite};
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode, Stdio};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

static CANCELLED: AtomicBool = AtomicBool::new(false);

// ─── Defaults ────────────────────────────────────────────────────────────────

const GENOME_DIR: &str =
    "/home/cml/humandb/transcriptomeindex/ensembl113/star_hg38_101bp_index";
const GTF_FILE: &str =
    "/home/cml/humandb/transcriptomeindex/ensembl113/Homo_sapiens.GRCh38.113.gtf";
const STAR_ENV: &str = "/home/cml/miniforge3/envs/star";
const RSEQC_ENV: &str = "/home/cml/miniforge3/envs/RSeQC";
const DEEPTOOLS_ENV: &str = "/home/cml/miniforge3/envs/deeptools";
const SAMTOOLS_BIN: &str = "samtools"; // resolved from $PATH by default

const REFRESH_INTERVAL: Duration = Duration::from_millis(100);

// ─── Help ────────────────────────────────────────────────────────────────────

fn usage() {
    eprintln!();
    let w = terminal::size().map(|(c, _)| c as usize).unwrap_or(80);
    let sep = "═".repeat(w);
    eprintln!("{sep}");
    eprintln!("{:^width$}", format!("STAR-RSeQC v{}", env!("CARGO_PKG_VERSION")), width = w);
    eprintln!("{:^width$}", "RNA-seq 2-pass alignment + quality control pipeline", width = w);
    eprintln!("{sep}");
    eprintln!();
    eprintln!("QUICK START:");
    eprintln!("    star-rseqc ./fastq/              # Run on current directory");
    eprintln!("    star-rseqc /data/paired -j 4     # Run with 4 parallel jobs");
    eprintln!();
    eprintln!("USAGE:");
    eprintln!("    star-rseqc <FASTQ_DIR> [OPTIONS]");
    eprintln!();
    eprintln!("DESCRIPTION:");
    eprintln!("    Complete RNA-seq 2-phase pipeline for paired-end reads:");
    eprintln!();
    eprintln!("    Phase 1: STAR 2-pass alignment with chimeric detection");
    eprintln!("       ✓ Produces coordinate-sorted BAM + transcriptome BAM");
    eprintln!("       ✓ Gene-level counts (ENCODE-compliant parameters)");
    eprintln!("       ✓ Chimeric junctions for fusion detection");
    eprintln!();
    eprintln!("    Phase 2: deeptools BAM → bigwig conversion");
    eprintln!("       ✓ bamCoverage           - Convert sorted BAM to bigwig (binSize 10)");
    eprintln!();
    eprintln!("    Phase 3: RSeQC quality control analysis (all 3 commands in parallel)");
    eprintln!("       ✓ infer_experiment.py   - Detect library strandedness");
    eprintln!("       ✓ read_distribution.py  - Genomic feature distribution analysis");
    eprintln!("       ✓ geneBody_coverage2.py - Check 5'-to-3' coverage uniformity (bigwig)");
    eprintln!("       ✓ ggplot2 PDF plots     - Curves + heatmap from coverage data");
    eprintln!();
    eprintln!("KEY FEATURES:");
    eprintln!("    • Three-phase pipeline: STAR → deeptools → RSeQC");
    eprintln!("    • Ultimate parallelism: all 3 RSeQC commands run simultaneously per sample");
    eprintln!("    • SHA256-based resume awareness: per-step checkpoint verification");
    eprintln!("    • Full-screen progress TUI with real-time updates");
    eprintln!("    • Auto-detects system resources (RAM, CPUs)");
    eprintln!("    • Auto-generates BED12 from GTF (cached for speed)");
    eprintln!("    • Graceful cancellation with Ctrl+C");
    eprintln!("    • JSON + TSV pipeline summaries");
    eprintln!();
    eprintln!("ARGUMENTS:");
    eprintln!("    <FASTQ_DIR>               Input directory with *_1P.fastq.gz files");
    eprintln!("                              (paired R2 must be *_2P.fastq.gz)");
    eprintln!();
    eprintln!("RESOURCE OPTIONS:");
    eprintln!("    -o, --output <DIR>        Output directory [default: star-rseqc-results]");
    eprintln!("    -j, --jobs <N>            Default parallel jobs [default: auto-detected from RAM]");
    eprintln!("    --star-jobs <N>           STAR phase jobs [default: same as --jobs]");
    eprintln!("    --deeptools-jobs <N>      deeptools phase jobs [default: same as --jobs]");
    eprintln!("    --rseqc-jobs <N>          RSeQC phase jobs [default: same as --jobs]");
    eprintln!("    -t, --threads <N>         Threads per job [default: auto-detected from CPUs]");
    eprintln!("    --bam-sort-ram <BYTES>    BAM sort RAM [default: auto-detected]");
    eprintln!();
    eprintln!("REFERENCE FILES:");
    eprintln!("    --genome-dir <DIR>        STAR genome index");
    eprintln!("                              [default: $STAR_RSEQC_GENOME_DIR or {}]", GENOME_DIR);
    eprintln!("    --gtf <FILE>              GTF annotation file");
    eprintln!("                              [default: $STAR_RSEQC_GTF or {}]", GTF_FILE);
    eprintln!("    --bed <FILE>              Pre-computed BED12 file (else auto-generated)");
    eprintln!();
    eprintln!("ENVIRONMENT & TOOLS:");
    eprintln!("    --star-env <DIR>          STAR environment path");
    eprintln!("                              [default: $STAR_RSEQC_STAR_ENV or {}]", STAR_ENV);
    eprintln!("    --rseqc-env <DIR>         RSeQC environment path");
    eprintln!("                              [default: $STAR_RSEQC_RSEQC_ENV or {}]", RSEQC_ENV);
    eprintln!("    --deeptools-env <DIR>     deeptools environment path (for bamCoverage)");
    eprintln!("                              [default: $STAR_RSEQC_DEEPTOOLS_ENV or {}]", DEEPTOOLS_ENV);
    eprintln!();
    eprintln!("ENVIRONMENT VARIABLE OVERRIDES:");
    eprintln!("    STAR_RSEQC_GENOME_DIR     Default for --genome-dir");
    eprintln!("    STAR_RSEQC_GTF            Default for --gtf");
    eprintln!("    STAR_RSEQC_STAR_ENV       Default for --star-env");
    eprintln!("    STAR_RSEQC_RSEQC_ENV      Default for --rseqc-env");
    eprintln!("    STAR_RSEQC_DEEPTOOLS_ENV  Default for --deeptools-env");
    eprintln!("    --samtools <PATH>         samtools binary path");
    eprintln!("                              [default: {}]", SAMTOOLS_BIN);
    eprintln!();
    eprintln!("WORKFLOW OPTIONS:");
    eprintln!("    --skip-alignment          QC only (skip STAR, run on existing BAMs)");
    eprintln!("    --dry-run                 Preview samples without running");
    eprintln!("    --clean-sam               Delete chimeric SAM files after successful run");
    eprintln!("    --operator <NAME>         Operator name for audit trail [default: $USER]");
    eprintln!();
    eprintln!("OTHER:");
    eprintln!("    -h, --help                Display this help message");
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
    eprintln!("    STAR index : {}", GENOME_DIR);
    eprintln!("    GTF        : {}", GTF_FILE);
    eprintln!();
    eprintln!("TOOL ENVIRONMENTS:");
    eprintln!("    STAR       : {}/bin/STAR", STAR_ENV);
    eprintln!("    samtools   : {}", SAMTOOLS_BIN);
    eprintln!("    RSeQC      : {}/bin/infer_experiment.py", RSEQC_ENV);
    eprintln!("    deeptools  : {}/bin/bamCoverage", DEEPTOOLS_ENV);
    eprintln!("    Rscript    : system PATH (optional, for ggplot2 PDF plots)");
    eprintln!();
    eprintln!("OUTPUT STRUCTURE:");
    eprintln!("    <output>/");
    eprintln!("      star/                          STAR alignment output per sample");
    eprintln!("        <sample>_Aligned.sortedByCoord.out.bam");
    eprintln!("        <sample>_Aligned.toTranscriptome.out.bam");
    eprintln!("        <sample>_ReadsPerGene.out.tab");
    eprintln!("        <sample>_Chimeric.out.junction");
    eprintln!("        <sample>_Log.final.out");
    eprintln!("      bigwig/                        deeptools bigwig files (Phase 2)");
    eprintln!("        <sample>.bw");
    eprintln!("      qc/                            RSeQC quality control output (Phase 3)");
    eprintln!("        <sample>.strand.txt");
    eprintln!("        <sample>.read_distribution.txt");
    eprintln!("        <sample>.geneBodyCoverage.txt       (coverage data)");
    eprintln!("        <sample>.geneBodyCoverage_plot.r    (RSeQC R script)");
    eprintln!("        <sample>.geneBodyCoverage.pdf       (RSeQC basic plot)");
    eprintln!("        <sample>.geneBodyCoverage.curves.pdf   (ggplot2 curves*)");
    eprintln!("        <sample>.geneBodyCoverage.heatMap.pdf  (ggplot2 heatmap*)");
    eprintln!("        * requires Rscript + ggplot2 (mandatory, included in checkpoint)");
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
    eprintln!();
    eprintln!("    PDF plots require Rscript + ggplot2 in PATH.");
    eprintln!("    Install: apt-get install r-base && Rscript -e 'install.packages(\"ggplot2\")'");
}

// ─── System resource detection ───────────────────────────────────────────────

fn read_available_ram() -> u64 {
    let Ok(file) = std::fs::File::open("/proc/meminfo") else { return 0; };
    let reader = BufReader::new(file);
    for line in reader.lines().flatten() {
        if let Some(rest) = line.strip_prefix("MemAvailable:") {
            let kib: u64 = rest.split_whitespace()
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
    let cpus = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
    (ram, cpus)
}

fn auto_config_resources(available_ram: u64, total_cpus: usize) -> (usize, usize, u64) {
    const PER_JOB_RAM: u64 = 38_000_000_000; // 32 GB STAR genome index + 6 GB BAM sort
    const BAM_SORT_MAX: u64 = 6_000_000_000; // ideal BAM sort allocation per job
    const OS_BUFFER: u64 = 2_000_000_000; // reserve for OS
    let usable = available_ram.saturating_sub(OS_BUFFER);
    let jobs = ((usable / PER_JOB_RAM) as usize).max(1);
    let threads = (total_cpus / jobs).max(1);
    // Cap BAM sort RAM at 25% of usable RAM divided across jobs, so we never
    // allocate more sort RAM than the system can provide on constrained machines.
    let bam_sort_ram = BAM_SORT_MAX.min(usable / (jobs as u64 * 4)).max(1_000_000_000);
    (jobs, threads, bam_sort_ram)
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
    deeptools_env: PathBuf,
    samtools: PathBuf,
    threads_per_sample: usize,
    parallel_jobs: usize,
    parallel_star_jobs: usize,
    parallel_rseqc_jobs: usize,
    parallel_deeptools_jobs: usize,
    bam_sort_ram: u64,
    skip_alignment: bool,
    dry_run: bool,
    clean_sam: bool,
    operator: String,
    resources_auto: bool,
}

/// Return value from parse_args:
///   Ok(Config)  — proceed normally
///   Err(true)   — --help was printed; exit 0
///   Err(false)  — argument error; exit 1
fn parse_args() -> Result<Config, bool> {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        usage();
        return Err(false);
    }

    let mut fastq_dir: Option<PathBuf> = None;
    let mut output_dir = PathBuf::from("star-rseqc-results");
    // Environment variables override compiled-in defaults, allowing deployment on
    // different machines without recompiling. Explicit flags override env vars.
    let mut genome_dir = PathBuf::from(
        env::var("STAR_RSEQC_GENOME_DIR").unwrap_or_else(|_| GENOME_DIR.to_string()));
    let mut gtf = PathBuf::from(
        env::var("STAR_RSEQC_GTF").unwrap_or_else(|_| GTF_FILE.to_string()));
    let mut bed: Option<PathBuf> = None;
    let mut star_env = PathBuf::from(
        env::var("STAR_RSEQC_STAR_ENV").unwrap_or_else(|_| STAR_ENV.to_string()));
    let mut rseqc_env = PathBuf::from(
        env::var("STAR_RSEQC_RSEQC_ENV").unwrap_or_else(|_| RSEQC_ENV.to_string()));
    let mut deeptools_env = PathBuf::from(
        env::var("STAR_RSEQC_DEEPTOOLS_ENV").unwrap_or_else(|_| DEEPTOOLS_ENV.to_string()));
    let mut samtools = PathBuf::from(SAMTOOLS_BIN);
    let mut threads_per_sample: Option<usize> = None;
    let mut parallel_jobs: Option<usize> = None;
    let mut parallel_star_jobs: Option<usize> = None;
    let mut parallel_rseqc_jobs: Option<usize> = None;
    let mut parallel_deeptools_jobs: Option<usize> = None;
    let mut bam_sort_ram: Option<u64> = None;
    let mut skip_alignment = false;
    let mut dry_run = false;
    let mut clean_sam = false;
    let mut operator: Option<String> = None;

    let mut i = 1;

    // Helper: consume the next argument as a value for a flag
    macro_rules! next_val {
        ($flag:expr) => {{
            i += 1;
            if i >= args.len() {
                eprintln!("Error: {} requires a value.", $flag);
                return Err(false);
            }
            &args[i]
        }};
    }
    // Helper: parse a usize flag value without calling process::exit
    macro_rules! parse_usize {
        ($flag:expr, $val:expr) => {
            match $val.parse::<usize>() {
                Ok(v) => v,
                Err(_) => {
                    eprintln!("Error: invalid value for {}: '{}' (expected a positive integer)", $flag, $val);
                    return Err(false);
                }
            }
        };
    }
    // Helper: parse a u64 flag value without calling process::exit
    macro_rules! parse_u64 {
        ($flag:expr, $val:expr) => {
            match $val.parse::<u64>() {
                Ok(v) => v,
                Err(_) => {
                    eprintln!("Error: invalid value for {}: '{}' (expected a positive integer)", $flag, $val);
                    return Err(false);
                }
            }
        };
    }

    while i < args.len() {
        match args[i].as_str() {
            "-h" | "--help" => {
                usage();
                return Err(true); // help printed; exit 0
            }
            "-o" | "--output" => {
                output_dir = PathBuf::from(next_val!("-o/--output"));
            }
            "-j" | "--jobs" => {
                let v = next_val!("-j/--jobs").clone();
                parallel_jobs = Some(parse_usize!("--jobs", v));
            }
            "--star-jobs" => {
                let v = next_val!("--star-jobs").clone();
                parallel_star_jobs = Some(parse_usize!("--star-jobs", v));
            }
            "--rseqc-jobs" => {
                let v = next_val!("--rseqc-jobs").clone();
                parallel_rseqc_jobs = Some(parse_usize!("--rseqc-jobs", v));
            }
            "--deeptools-jobs" => {
                let v = next_val!("--deeptools-jobs").clone();
                parallel_deeptools_jobs = Some(parse_usize!("--deeptools-jobs", v));
            }
            "-t" | "--threads" => {
                let v = next_val!("-t/--threads").clone();
                threads_per_sample = Some(parse_usize!("--threads", v));
            }
            "--genome-dir" => {
                genome_dir = PathBuf::from(next_val!("--genome-dir"));
            }
            "--gtf" => {
                gtf = PathBuf::from(next_val!("--gtf"));
            }
            "--bed" => {
                bed = Some(PathBuf::from(next_val!("--bed")));
            }
            "--star-env" => {
                star_env = PathBuf::from(next_val!("--star-env"));
            }
            "--rseqc-env" => {
                rseqc_env = PathBuf::from(next_val!("--rseqc-env"));
            }
            "--deeptools-env" => {
                deeptools_env = PathBuf::from(next_val!("--deeptools-env"));
            }
            "--samtools" => {
                samtools = PathBuf::from(next_val!("--samtools"));
            }
            "--bam-sort-ram" => {
                let v = next_val!("--bam-sort-ram").clone();
                bam_sort_ram = Some(parse_u64!("--bam-sort-ram", v));
            }
            "--skip-alignment" => skip_alignment = true,
            "--dry-run" => dry_run = true,
            "--clean-sam" => clean_sam = true,
            "--operator" => {
                operator = Some(next_val!("--operator").clone());
            }
            other => {
                if other.starts_with('-') {
                    eprintln!("Unknown option: {other}");
                    eprintln!("Run with -h for help.");
                    return Err(false);
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
            return Err(false);
        }
    };

    // Resolve resource parameters — auto-detect each missing value individually.
    // resources_auto is only true when everything was auto-detected (for display purposes).
    let resources_auto = parallel_jobs.is_none() && threads_per_sample.is_none() && bam_sort_ram.is_none();
    let (parallel_jobs, threads_per_sample, bam_sort_ram) = {
        let (avail_ram, total_cpus) = detect_system_resources();
        let (aj, at, ab) = auto_config_resources(avail_ram, total_cpus);
        (
            parallel_jobs.unwrap_or(aj),
            threads_per_sample.unwrap_or(at),
            bam_sort_ram.unwrap_or(ab),
        )
    };

    // Phase-specific job counts default to overall parallel_jobs if not specified
    let parallel_star_jobs = parallel_star_jobs.unwrap_or(parallel_jobs);
    let parallel_rseqc_jobs = parallel_rseqc_jobs.unwrap_or(parallel_jobs);
    let parallel_deeptools_jobs = parallel_deeptools_jobs.unwrap_or(parallel_jobs);

    // If --star-jobs overrides the auto-detected job count, recompute bam_sort_ram
    // so we don't hand STAR 6 GB sort RAM while running 8 parallel jobs on a 64 GB machine.
    let bam_sort_ram = if bam_sort_ram == {
        let (_, _, ab) = {
            let (avail_ram, total_cpus) = detect_system_resources();
            auto_config_resources(avail_ram, total_cpus)
        };
        ab
    } && parallel_star_jobs != parallel_jobs {
        // User overrode star-jobs but left bam_sort_ram on auto — recalculate for actual star parallelism
        let (avail_ram, _) = detect_system_resources();
        const BAM_SORT_MAX: u64 = 6_000_000_000;
        const OS_BUFFER: u64 = 2_000_000_000;
        let usable = avail_ram.saturating_sub(OS_BUFFER);
        BAM_SORT_MAX
            .min(usable / (parallel_star_jobs as u64 * 4))
            .max(1_000_000_000)
    } else {
        bam_sort_ram
    };

    if parallel_jobs == 0 {
        eprintln!("Error: --jobs must be >= 1");
        return Err(false);
    }
    if parallel_star_jobs == 0 {
        eprintln!("Error: --star-jobs must be >= 1");
        return Err(false);
    }
    if parallel_rseqc_jobs == 0 {
        eprintln!("Error: --rseqc-jobs must be >= 1");
        return Err(false);
    }
    if parallel_deeptools_jobs == 0 {
        eprintln!("Error: --deeptools-jobs must be >= 1");
        return Err(false);
    }
    if threads_per_sample == 0 {
        eprintln!("Error: --threads must be >= 1");
        return Err(false);
    }

    Ok(Config {
        fastq_dir,
        output_dir,
        genome_dir,
        gtf,
        bed,
        star_env,
        rseqc_env,
        deeptools_env,
        samtools,
        threads_per_sample,
        parallel_jobs,
        parallel_star_jobs,
        parallel_rseqc_jobs,
        parallel_deeptools_jobs,
        bam_sort_ram,
        skip_alignment,
        dry_run,
        clean_sam,
        operator: {
            let raw = operator.unwrap_or_else(|| {
                std::env::var("USER")
                    .or_else(|_| std::env::var("LOGNAME"))
                    .unwrap_or_else(|_| "unknown".to_string())
            });
            // Strip control characters and limit length so the audit JSON is always well-formed.
            let sanitized: String = raw.chars()
                .filter(|c| !c.is_control())
                .take(128)
                .collect();
            if sanitized.is_empty() { "unknown".to_string() } else { sanitized }
        },
        resources_auto,
    })
}

fn is_cancelled() -> bool {
    CANCELLED.load(Ordering::SeqCst)
}

// ─── Sample types ────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
struct Sample {
    name: String,
    r1: PathBuf,
    r2: PathBuf,
}

#[derive(Debug, Clone)]
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
    recent_events: Mutex<VecDeque<String>>,
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
            recent_events: Mutex::new(VecDeque::new()),
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
            events.push_back(msg);
            if events.len() > 100 {
                events.pop_front();
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
        hasher.update(path.file_name().unwrap_or_default().to_string_lossy().as_bytes());
        let size = fs::metadata(path).map(|m| m.len()).unwrap_or(0);
        if path.exists() {
            if size > 0 {
                if let Ok(bytes) = sha256_file(path) {
                    hasher.update(&bytes);
                }
            } else {
                // Zero-byte file: exists but empty (e.g. disk-full crash).
                // Use a distinct sentinel so it hashes differently from truly absent files.
                hasher.update(b"__ZERO_BYTES__");
            }
        } else {
            hasher.update(b"__MISSING__");
        }
    }
    format!("{:x}", hasher.finalize())
}

/// Per-step SHA256 digests — one field per individual tool output.
struct SampleDigests {
    star:     String,  // Phase 1: BAM + index + logs
    deeptools: String, // Phase 2: bigwig
    infer:    String,  // Phase 3a: strand.txt
    rdist:    String,  // Phase 3b: read_distribution.txt
    genebody: String,  // Phase 3c: geneBodyCoverage files + PDFs
}

/// Step identifiers for per-step checkpoint files.
const STEP_STAR: &str = "star";
const STEP_DEEPTOOLS: &str = "deeptools";
const STEP_INFER: &str = "infer";
const STEP_RDIST: &str = "rdist";
const STEP_GENEBODY: &str = "genebody";

/// Compute the SHA256 digest for a single step's output files.
/// All files are hashed unconditionally — missing or empty files are included
/// with sentinel values so the digest is always deterministic.
///
fn sha256_step(output_dir: &Path, sample_name: &str, step: &str) -> String {
    let star_dir = output_dir.join("star");
    let qc_dir   = output_dir.join("qc");
    let bw_dir   = output_dir.join("bigwig");

    match step {
        STEP_STAR => {
            // Hash log files + BAI + secondary outputs. Full BAM content hash is
            // intentionally omitted (multi-GB); the BAM's existence is captured via
            // the BAI and log entries which change whenever STAR reruns.
            sha256_file_list(&[
                star_dir.join(format!("{sample_name}_Log.out")),
                star_dir.join(format!("{sample_name}_Log.progress.out")),
                star_dir.join(format!("{sample_name}_Log.final.out")),
                star_dir.join(format!("{sample_name}_Aligned.sortedByCoord.out.bam.bai")),
                star_dir.join(format!("{sample_name}_Aligned.toTranscriptome.out.bam")),
                star_dir.join(format!("{sample_name}_ReadsPerGene.out.tab")),
                star_dir.join(format!("{sample_name}_SJ.out.tab")),
                star_dir.join(format!("{sample_name}_Chimeric.out.junction")),
                // _Chimeric.out.sam is intentionally excluded: --clean-sam may delete it
                // after alignment and its removal must not invalidate the STAR checkpoint.
            ])
        }
        STEP_DEEPTOOLS => {
            sha256_file_list(&[bw_dir.join(format!("{sample_name}.bw"))])
        }
        STEP_INFER => {
            sha256_file_list(&[qc_dir.join(format!("{sample_name}.strand.txt"))])
        }
        STEP_RDIST => {
            sha256_file_list(&[qc_dir.join(format!("{sample_name}.read_distribution.txt"))])
        }
        STEP_GENEBODY => {
            sha256_file_list(&[
                qc_dir.join(format!("{sample_name}.geneBodyCoverage.txt")),
                qc_dir.join(format!("{sample_name}.geneBodyCoverage_plot.r")),
                qc_dir.join(format!("{sample_name}.geneBodyCoverage.pdf")),
                qc_dir.join(format!("{sample_name}.geneBodyCoverage.curves.r")),
                qc_dir.join(format!("{sample_name}.geneBodyCoverage.curves.pdf")),
                qc_dir.join(format!("{sample_name}.geneBodyCoverage.heatmap.r")),
                qc_dir.join(format!("{sample_name}.geneBodyCoverage.heatMap.pdf")),
            ])
        }
        _ => "UNKNOWN_STEP".to_string(),
    }
}

/// Path for a per-step checkpoint file.
fn step_checkpoint_path(output_dir: &Path, sample_name: &str, step: &str) -> PathBuf {
    checkpoint_dir(output_dir).join(format!("{sample_name}.{step}.sha256"))
}

/// Write a per-step checkpoint file containing the SHA256 digest.
/// Called immediately after a step completes successfully.
fn write_step_checkpoint(output_dir: &Path, sample_name: &str, step: &str) -> Result<String, String> {
    let dir = checkpoint_dir(output_dir);
    fs::create_dir_all(&dir)
        .map_err(|e| format!("Failed to create checkpoint dir {}: {e}", dir.display()))?;
    let digest = sha256_step(output_dir, sample_name, step);
    let path = step_checkpoint_path(output_dir, sample_name, step);
    atomic_write(&path, digest.as_bytes())
        .map_err(|e| format!("Failed to write checkpoint {}: {e}", path.display()))?;
    Ok(digest)
}

/// Read a saved per-step checkpoint digest. Returns None if the file doesn't exist.
fn read_step_checkpoint(output_dir: &Path, sample_name: &str, step: &str) -> Option<String> {
    let path = step_checkpoint_path(output_dir, sample_name, step);
    fs::read_to_string(&path).ok().map(|s| s.trim().to_string())
}

/// Check whether a step is complete: the checkpoint file must exist AND its
/// stored digest must match the current output hash. This is the SOLE resume
/// authority — file existence/size is never checked directly.
fn is_step_done(output_dir: &Path, sample_name: &str, step: &str) -> bool {
    match read_step_checkpoint(output_dir, sample_name, step) {
        Some(saved) => {
            let current = sha256_step(output_dir, sample_name, step);
            saved == current
        }
        None => false,
    }
}

/// Delete a step's checkpoint file (used before rerunning a step).
fn remove_step_checkpoint(output_dir: &Path, sample_name: &str, step: &str) {
    let path = step_checkpoint_path(output_dir, sample_name, step);
    let _ = fs::remove_file(&path);
}

/// Write `data` to `path` atomically: write to a PID-tagged temp file in the same
/// directory, then rename. Prevents partial-write corruption on crash.
/// Using the same directory for src and dst keeps the rename within one filesystem,
/// which is required for POSIX rename(2) atomicity guarantees.
fn atomic_write(path: &Path, data: &[u8]) -> io::Result<()> {
    let parent = path.parent().unwrap_or(Path::new("."));
    let name = path.file_name().unwrap_or_default().to_string_lossy();
    // Include both PID and thread ID so concurrent threads writing different
    // files in the same directory never share a temp filename.
    let tid = {
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};
        let mut h = DefaultHasher::new();
        std::thread::current().id().hash(&mut h);
        h.finish()
    };
    let tmp = parent.join(format!(".{}.tmp{}-{:x}", name, std::process::id(), tid));
    {
        let mut f = File::create(&tmp)?;
        f.write_all(data)?;
        f.sync_all()?; // flush to disk before rename — crash-safe checkpoint
    }
    fs::rename(&tmp, path).map_err(|e| {
        let _ = fs::remove_file(&tmp);
        e
    })
}

/// Check resume status for a sample's 3-phase outputs using per-step checkpoint files.
/// Each step has its own .sha256 file written immediately after successful completion.
/// A step is "done" only if its checkpoint file exists AND the stored digest matches
/// the current output hash. No file existence/size checks — purely checkpoint-driven.
///
/// Returns:
///   AllDone       — all steps complete and intact → skip
///   Phase1Changed — STAR checkpoint missing/mismatched → redo all
///   Phase2Changed — deeptools checkpoint missing/mismatched → redo Phase 2 & 3
///   Phase3Changed — any RSeQC checkpoint missing/mismatched → redo Phase 3 only
///   NotDone       — no checkpoints found → process all phases
#[derive(Clone)]
enum ResumeStatus {
    AllDone,        // Phase 1 (STAR) ✓, Phase 2 (deeptools) ✓, Phase 3 (RSeQC) ✓
    Phase1Changed,  // STAR checkpoint missing/mismatched → redo all
    Phase2Changed,  // deeptools checkpoint missing/mismatched → redo Phase 2 & 3
    Phase3Changed,  // RSeQC checkpoint(s) missing/mismatched → redo Phase 3 only
    NotDone,        // No checkpoints found
}

fn check_resume(output_dir: &Path, sample_name: &str) -> ResumeStatus {
    let star_ok      = is_step_done(output_dir, sample_name, STEP_STAR);
    let deeptools_ok = is_step_done(output_dir, sample_name, STEP_DEEPTOOLS);
    let infer_ok     = is_step_done(output_dir, sample_name, STEP_INFER);
    let rdist_ok     = is_step_done(output_dir, sample_name, STEP_RDIST);
    let genebody_ok  = is_step_done(output_dir, sample_name, STEP_GENEBODY);
    let phase3_ok    = infer_ok && rdist_ok && genebody_ok;

    // If none of the checkpoints exist at all, treat as NotDone
    let any_exists = star_ok || deeptools_ok || infer_ok || rdist_ok || genebody_ok;
    if !any_exists {
        // Also check if any checkpoint files exist but are mismatched (legacy migration)
        let has_any_file = [STEP_STAR, STEP_DEEPTOOLS, STEP_INFER, STEP_RDIST, STEP_GENEBODY]
            .iter()
            .any(|step| step_checkpoint_path(output_dir, sample_name, step).exists());
        if !has_any_file {
            // Also try to migrate from legacy monolithic checkpoint format
            let legacy = checkpoint_dir(output_dir).join(format!("{sample_name}.sha256"));
            if legacy.exists() {
                // Legacy checkpoint found — remove it and reprocess from scratch
                let _ = fs::remove_file(&legacy);
            }
            return ResumeStatus::NotDone;
        }
    }

    match (star_ok, deeptools_ok, phase3_ok) {
        (true, true, true) => ResumeStatus::AllDone,
        (true, true, false) => ResumeStatus::Phase3Changed,
        (true, false, _) => ResumeStatus::Phase2Changed,
        (false, _, _) => ResumeStatus::Phase1Changed,
    }
}

/// Number of threads to use for parallel resume checks.
/// Defined once here so the display log and the actual spawn count always agree.
fn resume_check_threads(n: usize) -> usize { n.min(32).max(1) }

/// Parallel resume detection: check all samples simultaneously using thread pool.
/// Samples are shared via Arc to avoid cloning the full list for each thread.
fn check_resume_all_parallel(output_dir: &Path, samples: &[Sample]) -> Vec<(String, ResumeStatus)> {
    let results: Arc<Mutex<Vec<(String, ResumeStatus)>>> = Arc::new(Mutex::new(Vec::new()));
    let next_idx = Arc::new(AtomicUsize::new(0));
    let num_threads = resume_check_threads(samples.len());
    // Share samples and output_dir across threads without cloning each per thread.
    let shared_samples: Arc<Vec<Sample>> = Arc::new(samples.to_vec());
    let shared_od: Arc<PathBuf> = Arc::new(output_dir.to_path_buf());

    let mut handles = vec![];

    for _ in 0..num_threads {
        let next_arc = Arc::clone(&next_idx);
        let res_arc = Arc::clone(&results);
        let samps = Arc::clone(&shared_samples);
        let od = Arc::clone(&shared_od);

        let handle = std::thread::spawn(move || {
            loop {
                let idx = next_arc.fetch_add(1, Ordering::Relaxed);
                if idx >= samps.len() {
                    break;
                }
                let sample = &samps[idx];
                let status = check_resume(&od, &sample.name);
                res_arc.lock().unwrap_or_else(|e| e.into_inner()).push((sample.name.clone(), status));
            }
        });

        handles.push(handle);
    }

    for handle in handles {
        if handle.join().is_err() {
            eprintln!("Warning: SHA256 worker thread panicked — affected samples will be reprocessed");
        }
    }

    let final_results = results.lock().unwrap_or_else(|e| e.into_inner()).clone();
    final_results
}

// ─── Sample discovery ────────────────────────────────────────────────────────

fn discover_samples(fastq_dir: &Path) -> Vec<Sample> {
    let pattern = fastq_dir
        .join("*_1P.fastq.gz")
        .to_string_lossy()
        .to_string();

    let mut samples = Vec::new();
    let mut seen = HashMap::new();

    const MAX_SAMPLES: usize = 5_000;
    let entries: Vec<_> = match glob(&pattern) {
        Ok(paths) => {
            let v: Vec<_> = paths.filter_map(|e| e.ok()).collect();
            if v.len() > MAX_SAMPLES {
                eprintln!("Warning: {} FASTQ R1 files found — this is unusually large. \
                    Check that you specified the correct FASTQ directory (max {} samples supported).",
                    v.len(), MAX_SAMPLES);
            }
            v
        }
        Err(_) => return samples,
    };

    for r1 in entries {
        let r1_name = r1.file_name().unwrap().to_string_lossy().to_string();

        let sample_name = match r1_name.strip_suffix("_1P.fastq.gz") {
            Some(n) => n.to_string(),
            None => continue,
        };

        // Reject names that could traverse paths or inject shell characters.
        // Only alphanumerics, underscore, hyphen, and dot are permitted.
        if !sample_name.chars().all(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '-' | '.')) {
            eprintln!("Warning: skipping {} — sample name contains characters outside [A-Za-z0-9_.-]", sample_name);
            continue;
        }
        if sample_name.is_empty() || sample_name.starts_with('.') {
            eprintln!("Warning: skipping sample with empty or dot-prefixed name ({})", r1_name);
            continue;
        }

        let r2_name = format!("{}_2P.fastq.gz", sample_name);
        let r2 = r1.parent().unwrap().join(&r2_name);

        if !r2.exists() {
            eprintln!("Warning: skipping {} — R2 not found ({})", sample_name, r2.display());
            continue;
        }
        let r1_size = fs::metadata(&r1).map(|m| m.len()).unwrap_or(0);
        let r2_size = fs::metadata(&r2).map(|m| m.len()).unwrap_or(0);
        if r1_size == 0 || r2_size == 0 {
            eprintln!("Warning: skipping {} — FASTQ file is empty (R1={}B R2={}B)", sample_name, r1_size, r2_size);
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
        let start: u64 = match fields[3].parse::<u64>() {
            Ok(v) => v.saturating_sub(1),
            Err(_) => continue, // malformed coordinate — skip exon
        };
        let end: u64 = match fields[4].parse::<u64>() {
            Ok(v) => v,
            Err(_) => continue, // malformed coordinate — skip exon
        };
        if end <= start { continue; } // zero-length or inverted exon — skip
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

    // Write to a temp file first; rename atomically on success so a kill mid-write
    // never leaves a partial BED that gets reused on the next run.
    // If the process is OOM-killed while transcripts is in memory, this tmp file
    // will be left on disk. We remove any pre-existing tmp before starting so that
    // a stale tmp from a previous OOM crash does not accumulate.
    let tmp_path = bed_path.with_extension("bed12.tmp");
    let _ = fs::remove_file(&tmp_path); // remove stale tmp if present (ignore error if absent)
    let out_file = File::create(&tmp_path)
        .map_err(|e| format!("Cannot create BED tmp {}: {}", tmp_path.display(), e))?;
    let mut writer = BufWriter::new(out_file);

    // Sort by chromosome then start position for deterministic output
    let mut tx_vec: Vec<(String, (String, String, Vec<(u64, u64)>))> = transcripts.into_iter().collect();
    tx_vec.sort_by(|a, b| {
        let (chrom_a, _, exons_a) = &a.1;
        let (chrom_b, _, exons_b) = &b.1;
        let start_a = exons_a.iter().map(|e| e.0).min().unwrap_or(0);
        let start_b = exons_b.iter().map(|e| e.0).min().unwrap_or(0);
        chrom_a.cmp(chrom_b).then(start_a.cmp(&start_b))
    });

    let mut count = 0usize;
    for (tx_id, (chrom, strand, ref mut exons)) in &mut tx_vec {
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
        let _ = fs::remove_file(&tmp_path);
        return Err("GTF→BED12 produced zero transcripts".to_string());
    }
    // Flush before rename to ensure all bytes are on disk
    writer.flush().map_err(|e| format!("BED flush error: {e}"))?;
    drop(writer);
    fs::rename(&tmp_path, bed_path)
        .map_err(|e| format!("Cannot rename BED tmp to {}: {e}", bed_path.display()))?;
    Ok(count)
}

// ─── Run command with cancellation ───────────────────────────────────────────

/// Send SIGKILL and wait up to 5 s for the process to exit.
/// After 5 s, stop waiting — the process may be in uninterruptible I/O sleep
/// (stalled NFS), in which case it will be reaped by the OS on process exit.
fn kill_with_timeout(child: &mut std::process::Child) {
    let _ = child.kill();
    let deadline = Instant::now() + Duration::from_secs(5);
    loop {
        match child.try_wait() {
            Ok(Some(_)) | Err(_) => return,
            Ok(None) if Instant::now() >= deadline => {
                eprintln!("Warning: child process did not exit within 5 s of SIGKILL (possible stalled I/O)");
                return;
            }
            Ok(None) => std::thread::sleep(Duration::from_millis(100)),
        }
    }
}

fn run_cancellable(mut cmd: Command) -> Result<bool, String> {
    let mut child = cmd.spawn().map_err(|e| format!("Failed to launch: {e}"))?;

    loop {
        match child.try_wait() {
            Ok(Some(status)) => return Ok(status.success()),
            Ok(None) => {
                if is_cancelled() {
                    kill_with_timeout(&mut child);
                    return Err("Cancelled".to_string());
                }
                std::thread::sleep(Duration::from_millis(200));
            }
            Err(e) => return Err(format!("Wait error: {e}")),
        }
    }
}

/// Shared implementation for capture variants.
/// `combine` = true → stdout + stderr concatenated (for infer_experiment.py which writes to stderr).
/// `combine` = false → stdout only, stderr drained silently.
fn run_cancellable_capture_impl(mut cmd: Command, combine: bool) -> Result<Option<Vec<u8>>, String> {
    cmd.stdout(Stdio::piped()).stderr(Stdio::piped());
    let mut child = cmd.spawn().map_err(|e| format!("Failed to launch: {e}"))?;
    // Drain stdout AND stderr in background threads to prevent deadlock when
    // either pipe buffer fills (> ~64 KB on Linux). Both must be drained concurrently.
    let stdout_handle = child.stdout.take().map(|mut out| {
        std::thread::spawn(move || {
            let mut buf = Vec::new();
            let _ = out.read_to_end(&mut buf);
            buf
        })
    });
    let stderr_handle = child.stderr.take().map(|mut err| {
        std::thread::spawn(move || {
            let mut buf = Vec::new();
            let _ = err.read_to_end(&mut buf);
            buf
        })
    });
    loop {
        match child.try_wait() {
            Ok(Some(status)) => {
                let mut out_bytes = stdout_handle.and_then(|h| h.join().ok()).unwrap_or_default();
                let err_bytes = stderr_handle.and_then(|h| h.join().ok()).unwrap_or_default();
                if combine { out_bytes.extend_from_slice(&err_bytes); }
                return if status.success() { Ok(Some(out_bytes)) } else { Ok(None) };
            }
            Ok(None) => {
                if is_cancelled() {
                    kill_with_timeout(&mut child);
                    // Drain threads may block if the child is in uninterruptible sleep (stalled
                    // NFS). Join them with a timeout: spawn a joiner thread and give it 3 s.
                    let join_with_timeout = |h: std::thread::JoinHandle<Vec<u8>>| {
                        let (tx, rx) = std::sync::mpsc::channel();
                        std::thread::spawn(move || { let _ = tx.send(h.join()); });
                        let _ = rx.recv_timeout(Duration::from_secs(3));
                    };
                    if let Some(h) = stdout_handle { join_with_timeout(h); }
                    if let Some(h) = stderr_handle { join_with_timeout(h); }
                    return Err("Cancelled".to_string());
                }
                std::thread::sleep(Duration::from_millis(200));
            }
            Err(e) => {
                let _ = stdout_handle.and_then(|h| h.join().ok());
                let _ = stderr_handle.and_then(|h| h.join().ok());
                return Err(format!("Wait error: {e}"));
            }
        }
    }
}

/// Run a command with cancellation support, capturing stdout.
/// Returns Ok(Some(bytes)) on success, Ok(None) on non-zero exit, Err on cancel/spawn failure.
fn run_cancellable_capture(cmd: Command) -> Result<Option<Vec<u8>>, String> {
    run_cancellable_capture_impl(cmd, false)
}

/// Like run_cancellable_capture but returns stdout+stderr concatenated.
/// Use for tools like infer_experiment.py that write results to stderr.
fn run_cancellable_capture_combined(cmd: Command) -> Result<Option<Vec<u8>>, String> {
    run_cancellable_capture_impl(cmd, true)
}

/// Run `bin flag` and return the first non-empty output line, with a 5-second timeout.
/// Returns "unknown" if the command fails, times out, or produces no output.
/// The timeout prevents a broken conda environment from hanging the audit write at run end.
fn get_tool_version(bin: &Path, flag: &str) -> String {
    let bin = bin.to_path_buf();
    let flag = flag.to_string();
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let result = Command::new(&bin)
            .arg(&flag)
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .output()
            .ok()
            .map(|out| {
                let combined = [out.stdout.as_slice(), out.stderr.as_slice()].concat();
                String::from_utf8_lossy(&combined)
                    .lines()
                    .find(|l| !l.trim().is_empty())
                    .unwrap_or("unknown")
                    .trim()
                    .to_string()
            })
            .unwrap_or_else(|| "unknown".to_string());
        let _ = tx.send(result);
    });
    rx.recv_timeout(Duration::from_secs(5)).unwrap_or_else(|_| "unknown(timeout)".to_string())
}

fn make_log_stdio(log_dir: &Path, name: &str) -> Result<(Stdio, Stdio), String> {
    let log_path = log_dir.join(format!("{name}.log"));
    // Append mode: preserves the previous run's log on re-run so failure evidence is not lost.
    match std::fs::OpenOptions::new().create(true).append(true).open(&log_path) {
        Ok(f) => {
            let f2 = f.try_clone()
                .map_err(|e| format!("Failed to clone log file handle: {e}"))?;
            Ok((Stdio::from(f), Stdio::from(f2)))
        }
        Err(e) => Err(format!("Cannot open log file {}: {e}", log_path.display())),
    }
}

/// Generate R script for gene body coverage curves plot using ggplot2
fn create_curves_r_script(output_prefix: &str) -> Result<(), String> {
    let txt_file = format!("{}.geneBodyCoverage.txt", output_prefix);
    let pdf_file = format!("{}.geneBodyCoverage.curves.pdf", output_prefix);

    // Validate input file exists
    if !Path::new(&txt_file).exists() {
        return Err(format!("Input file not found: {}", txt_file));
    }

    let script = format!(
        r#"tryCatch({{
  # Load required packages
  if (!require('ggplot2', quietly=TRUE)) {{
    stop('ggplot2 package required for visualization. Install with: install.packages("ggplot2")')
  }}

  # Read data
  data <- read.table('{input}', header=TRUE, sep='\t')

  # Validate data
  if (nrow(data) == 0) {{
    stop('Input file is empty')
  }}
  if (!all(c('percentile', 'count') %in% names(data))) {{
    stop('Input file must contain percentile and count columns')
  }}

  # Create plot
  p <- ggplot(data, aes(x=percentile, y=count)) +
    geom_point(size=2, color='steelblue') +
    geom_line(color='steelblue', linewidth=0.8) +
    labs(
      title='Gene Body Coverage Distribution',
      x='Gene Body Percentile (5\' to 3\')',
      y='Average Coverage',
      subtitle='Coverage across gene body regions'
    ) +
    theme_bw() +
    theme(
      plot.title=element_text(size=14, face='bold'),
      plot.subtitle=element_text(size=11, color='grey40'),
      axis.title=element_text(size=12),
      panel.grid.minor=element_blank()
    )

  # Save plot
  ggsave('{output}', p, width=10, height=6, dpi=150)
}}, error = function(e) {{
  cat('Error generating curves plot:', conditionMessage(e), '\n', file=stderr())
  quit(status=1)
}})
"#,
        input = txt_file,
        output = pdf_file
    );

    let script_file = format!("{}.geneBodyCoverage.curves.r", output_prefix);
    std::fs::write(&script_file, script)
        .map_err(|e| format!("Failed to write curves R script: {e}"))
}

/// Generate R script for gene body coverage heatmap plot using ggplot2
fn create_heatmap_r_script(output_prefix: &str) -> Result<(), String> {
    let txt_file = format!("{}.geneBodyCoverage.txt", output_prefix);
    let pdf_file = format!("{}.geneBodyCoverage.heatMap.pdf", output_prefix);

    // Validate input file exists
    if !Path::new(&txt_file).exists() {
        return Err(format!("Input file not found: {}", txt_file));
    }

    let script = format!(
        r#"tryCatch({{
  # Load required packages
  if (!require('ggplot2', quietly=TRUE)) {{
    stop('ggplot2 package required for visualization. Install with: install.packages("ggplot2")')
  }}

  # Read data
  data <- read.table('{input}', header=TRUE, sep='\t')

  # Validate data
  if (nrow(data) == 0) {{
    stop('Input file is empty')
  }}
  if (!all(c('percentile', 'count') %in% names(data))) {{
    stop('Input file must contain percentile and count columns')
  }}

  # Normalize counts for heatmap visualization (0-1 scale)
  data$norm_count <- (data$count - min(data$count)) / (max(data$count) - min(data$count))

  # Create heatmap
  p <- ggplot(data, aes(x=percentile, y=1, fill=norm_count)) +
    geom_tile(height=0.8) +
    scale_fill_gradient(low='white', high='darkblue', name='Normalized\nCoverage') +
    scale_x_continuous(breaks=seq(0, 100, 10)) +
    labs(
      title='Gene Body Coverage Heatmap',
      x='Gene Body Percentile (5\' to 3\')',
      y='',
      subtitle='Intensity represents relative coverage level'
    ) +
    theme_bw() +
    theme(
      plot.title=element_text(size=14, face='bold'),
      plot.subtitle=element_text(size=11, color='grey40'),
      axis.title=element_text(size=12),
      axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      panel.grid=element_blank()
    )

  # Save plot
  ggsave('{output}', p, width=12, height=3, dpi=150)
}}, error = function(e) {{
  cat('Error generating heatmap plot:', conditionMessage(e), '\n', file=stderr())
  quit(status=1)
}})
"#,
        input = txt_file,
        output = pdf_file
    );

    let script_file = format!("{}.geneBodyCoverage.heatmap.r", output_prefix);
    std::fs::write(&script_file, script)
        .map_err(|e| format!("Failed to write heatmap R script: {e}"))
}

/// Run an R script via Rscript with proper error handling
fn run_r_script(rscript_path: &str) -> Result<(), String> {
    // Validate script exists
    if !Path::new(rscript_path).exists() {
        return Err(format!("R script not found: {}", rscript_path));
    }

    let mut cmd = Command::new("Rscript");
    cmd.arg("--vanilla").arg(rscript_path);

    match run_cancellable(cmd) {
        Ok(true) => Ok(()),
        Ok(false) => Err("Rscript exited with non-zero status".to_string()),
        Err(e) if e == "Cancelled" => Err("Cancelled".to_string()),
        Err(e) => Err(format!("Rscript error: {e}")),
    }
}

/// Validate that a PDF was successfully created
fn validate_pdf_created(pdf_path: &str) -> Result<(), String> {
    let path = Path::new(pdf_path);
    if !path.exists() {
        return Err(format!("PDF was not created: {}", pdf_path));
    }

    // Check that file is not empty (should be at least a few hundred bytes for PDF)
    match std::fs::metadata(path) {
        Ok(metadata) if metadata.len() > 100 => Ok(()),
        Ok(_) => Err(format!("PDF file appears empty or corrupted: {}", pdf_path)),
        Err(e) => Err(format!("Cannot read PDF file: {}: {}", pdf_path, e)),
    }
}

// ─── Pipeline steps ──────────────────────────────────────────────────────────

fn run_star_sample(
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
    let log_dir = config.output_dir.join("logs");

    let bam_path = star_dir.join(format!(
        "{}_Aligned.sortedByCoord.out.bam",
        sample.name
    ));
    let out_prefix = star_dir.join(format!("{}_", sample.name))
        .to_str()
        .ok_or_else(|| format!("{}: output path is not valid UTF-8", sample.name))?
        .to_string();
    let job_start = Instant::now();

    // ── Step 1: STAR Alignment ──
    if !config.skip_alignment {
        state.set_active(slot, &sample.name, "STAR alignment");

        // Always run STAR when this function is called — the SHA256 checkpoint is
        // the sole resume authority. Any pre-existing partial BAM is cleaned up first.
        {
            if bam_path.exists() {
                let bam_size = fs::metadata(&bam_path).map(|m| m.len()).unwrap_or(0);
                state.add_event(format!("  WARN  {} — removing existing BAM ({}B) before re-alignment", sample.name, bam_size));
                cleanup_partial_star(&star_dir, &sample.name);
            }
            let star_bin = config.star_env.join("bin/STAR");
            let (stdout_cfg, stderr_cfg) =
                make_log_stdio(&log_dir, &format!("{}.star", sample.name))?;

            let genome_dir_str = config.genome_dir.to_str()
                .ok_or_else(|| format!("{}: genome-dir path is not valid UTF-8", sample.name))?;
            let r1_str = sample.r1.to_str()
                .ok_or_else(|| format!("{}: R1 path is not valid UTF-8", sample.name))?;
            let r2_str = sample.r2.to_str()
                .ok_or_else(|| format!("{}: R2 path is not valid UTF-8", sample.name))?;
            let gtf_str = config.gtf.to_str()
                .ok_or_else(|| format!("{}: GTF path is not valid UTF-8", sample.name))?;

            let mut cmd = Command::new(&star_bin);
            cmd.args([
                "--runThreadN",
                &config.threads_per_sample.to_string(),
                "--genomeDir",
                genome_dir_str,
                "--readFilesIn",
                r1_str,
                r2_str,
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
                gtf_str,
            ])
            .stdout(stdout_cfg)
            .stderr(stderr_cfg);

            match run_cancellable(cmd) {
                Ok(true) => {
                    state.add_event(format!("  DONE  {} — STAR alignment", sample.name));
                }
                Ok(false) => {
                    cleanup_partial_star(&star_dir, &sample.name);
                    state.add_event(format!("  FAIL  {} — STAR alignment", sample.name));
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
    if !config.skip_alignment && !bam_path.exists() {
        state.add_event(format!("  FAIL  {} — BAM not found after alignment", sample.name));
        return Err(format!("{}: BAM not found: {}", sample.name, bam_path.display()));
    }

    // ── Step 2: samtools index ──
    // Under --skip-alignment the BAM comes from the user's existing data.
    // If no BAM is present, skip indexing silently (downstream phases will
    // filter this sample out via the BAM existence check in Phase 2).
    if !bam_path.exists() {
        state.add_event(format!("  SKIP  {} — no BAM present, skipping samtools index", sample.name));
        return Ok(());
    }
    state.update_step(slot, "samtools index");
    let bai_path = PathBuf::from(format!("{}.bai", bam_path.display()));
    // Treat a zero-byte BAI as absent — it means a previous index attempt crashed
    let bai_valid = bai_path.exists()
        && fs::metadata(&bai_path).map(|m| m.len() > 0).unwrap_or(false);
    if !bai_valid {
        let bam_str = bam_path.to_str()
            .ok_or_else(|| format!("{}: BAM path is not valid UTF-8", sample.name))?;
        let (stdout_cfg, stderr_cfg) =
            make_log_stdio(&log_dir, &format!("{}.samtools", sample.name))?;
        let mut cmd = Command::new(&config.samtools);
        cmd.args(["index", "-@", &config.threads_per_sample.to_string(), bam_str]);
        cmd.stdout(stdout_cfg).stderr(stderr_cfg);

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

    // ── Step 3: chr-prefix sanity check ──────────────────────────────────────
    // Warn if BAM chromosome names use "chr" prefix but BED12 does not (or vice
    // versa). This mismatch silently causes zero-coverage in all RSeQC tools.
    // Uses the resolved bed_path (works for both --bed and auto-generated BED12).
    if bed_path.exists() {
        let bam_hdr = Command::new(&config.samtools)
            .args(["view", "-H"])
            .arg(&bam_path)
            .output();
        if let Ok(hdr) = bam_hdr {
            let hdr_str = String::from_utf8_lossy(&hdr.stdout);
            let bam_has_chr = hdr_str.lines().any(|l| {
                l.starts_with("@SQ") && l.contains("\tSN:chr")
            });
            // Sample just the first data line in BED12
            let bed_has_chr = File::open(bed_path).ok()
                .and_then(|f| {
                    BufReader::new(f).lines().find_map(|l| {
                        l.ok().filter(|s| !s.is_empty() && !s.starts_with('#'))
                    })
                })
                .map(|l| l.starts_with("chr"))
                .unwrap_or(false);
            if bam_has_chr != bed_has_chr {
                state.add_event(format!(
                    "  WARN  {} — chr-prefix mismatch: BAM {} 'chr' prefix but BED12 {}; \
                     RSeQC results may be empty",
                    sample.name,
                    if bam_has_chr { "has" } else { "lacks" },
                    if bed_has_chr { "has it" } else { "lacks it" },
                ));
            }
        }
    }

    let dur = job_start.elapsed().as_secs_f64();
    state.record_duration(dur);
    Ok(())
}

/// Phase 2: deeptools — BAM → bigwig via bamCoverage
fn run_deeptools_phase(
    sample: &Sample,
    config: &Config,
    state: &ProgressState,
    slot: usize,
) -> Result<(), String> {
    if is_cancelled() {
        return Err("Cancelled".to_string());
    }

    let star_dir = config.output_dir.join("star");
    let bam_path = star_dir.join(format!("{}_Aligned.sortedByCoord.out.bam", sample.name));
    let bw_path = config.output_dir.join("bigwig").join(format!("{}.bw", sample.name));

    // bamCoverage requires a BAM index (.bai). Assert it exists with a clear diagnostic
    // so --skip-alignment users get an actionable error rather than a cryptic bamCoverage failure.
    let bai_path = PathBuf::from(format!("{}.bai", bam_path.display()));
    let bai_valid = bai_path.exists()
        && fs::metadata(&bai_path).map(|m| m.len() > 0).unwrap_or(false);
    if !bai_valid {
        return Err(format!(
            "{}: BAM index ({}) is missing or empty — bamCoverage requires a .bai index.\n\
             Run samtools index on the BAM first, or re-run without --skip-alignment.",
            sample.name, bai_path.display()
        ));
    }

    // Remove any pre-existing bigwig before running — SHA256 checkpoint is the sole
    // resume authority; a leftover bigwig from a prior run must not short-circuit this.
    if bw_path.exists() {
        let _ = fs::remove_file(&bw_path);
    }

    let t_start = Instant::now();

    state.set_active(slot, &sample.name, "bamCoverage (BAM -> bigwig)");
    let log_dir = config.output_dir.join("logs");
    let bam_coverage = config.deeptools_env.join("bin/bamCoverage");
    let bam_str = bam_path.to_str()
        .ok_or_else(|| format!("{}: BAM path is not valid UTF-8", sample.name))?;
    let bw_str = bw_path.to_str()
        .ok_or_else(|| format!("{}: bigwig path is not valid UTF-8", sample.name))?;
    let (stdout_cfg, stderr_cfg) =
        make_log_stdio(&log_dir, &format!("{}.bamcoverage", sample.name))?;
    let mut cmd = Command::new(&bam_coverage);
    cmd.args([
        "-b", bam_str,
        "-o", bw_str,
        "-p", &config.threads_per_sample.to_string(),
        "--binSize", "10",
    ]);
    cmd.stdout(stdout_cfg).stderr(stderr_cfg);

    match run_cancellable(cmd) {
        Ok(true) => {
            state.add_event(format!("  DONE  {} — bamCoverage", sample.name));
            state.record_duration(t_start.elapsed().as_secs_f64());
            Ok(())
        }
        Ok(false) => {
            let _ = fs::remove_file(&bw_path);
            state.add_event(format!("  FAIL  {} — bamCoverage exited non-zero", sample.name));
            Err(format!("{}: bamCoverage exited with error", sample.name))
        }
        Err(e) if e == "Cancelled" => {
            let _ = fs::remove_file(&bw_path);
            Err(e)
        }
        Err(e) => {
            let _ = fs::remove_file(&bw_path);
            state.add_event(format!("  FAIL  {} — bamCoverage: {e}", sample.name));
            Err(e)
        }
    }
}

/// Phase 3: RSeQC — infer_experiment + read_distribution + geneBody_coverage2 (all in parallel)
/// Requires bigwig from Phase 2 (bamCoverage) to be present.
fn run_rseqc_phase3(
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
    let bam_path = star_dir.join(format!("{}_Aligned.sortedByCoord.out.bam", sample.name));
    let bw_path = config.output_dir.join("bigwig").join(format!("{}.bw", sample.name));

    // SHA256 checkpoint is the sole resume authority — clean any pre-existing Phase 3
    // outputs so per-tool file_ok checks inside the threads start from a clean state.
    let strand_out = qc_dir.join(format!("{}.strand.txt", sample.name));
    let rdist_out = qc_dir.join(format!("{}.read_distribution.txt", sample.name));
    {
        let rseqc_python = config.rseqc_env.join("bin/python");
        state.set_active(slot, &sample.name, "infer + read_dist + geneBody_coverage2");

        let t_start = Instant::now();

        let infer_failed = Arc::new(AtomicBool::new(false));
        let rdist_failed = Arc::new(AtomicBool::new(false));
        let genebody_failed = Arc::new(AtomicBool::new(false));

        // Pre-compute UTF-8 path strings before entering thread::scope
        let bam_str = bam_path.to_str()
            .ok_or_else(|| format!("{}: BAM path is not valid UTF-8", sample.name))?
            .to_string();
        let bed_str = bed_path.to_str()
            .ok_or_else(|| format!("{}: BED path is not valid UTF-8", sample.name))?
            .to_string();
        let bw_str = bw_path.to_str()
            .ok_or_else(|| format!("{}: bigwig path is not valid UTF-8", sample.name))?
            .to_string();
        let genebody_prefix = qc_dir.join(&sample.name)
            .to_str()
            .ok_or_else(|| format!("{}: QC prefix path is not valid UTF-8", sample.name))?
            .to_string();

        // Run all 3 RSeQC steps in parallel via thread::scope
        {
            let infer_failed = Arc::clone(&infer_failed);
            let rdist_failed = Arc::clone(&rdist_failed);
            let genebody_failed = Arc::clone(&genebody_failed);
            std::thread::scope(|s| {
                // infer_experiment.py
                let t_infer = s.spawn(|| {
                    let script = config.rseqc_env.join("bin/infer_experiment.py");
                        let script_str = match script.to_str() {
                            Some(s) => s.to_string(),
                            None => {
                                infer_failed.store(true, Ordering::SeqCst);
                                state.add_event(format!("  FAIL  {} — infer_experiment script path is not valid UTF-8", sample.name));
                                return;
                            }
                        };
                        let mut cmd = Command::new(&rseqc_python);
                        cmd.args([
                            script_str.as_str(),
                            "-i", bam_str.as_str(),
                            "-r", bed_str.as_str(),
                        ]);
                        // infer_experiment.py writes results to stderr; capture both streams.
                        match run_cancellable_capture_combined(cmd) {
                            Ok(Some(combined)) => {
                                if let Err(e) = atomic_write(&strand_out, &combined) {
                                    state.add_event(format!("  FAIL  {} — infer_experiment write failed: {e}", sample.name));
                                    infer_failed.store(true, Ordering::SeqCst);
                                } else {
                                    state.add_event(format!("  DONE  {} — infer_experiment", sample.name));
                                }
                            }
                            Ok(None) => {
                                infer_failed.store(true, Ordering::SeqCst);
                                state.add_event(format!("  FAIL  {} — infer_experiment exited non-zero", sample.name));
                            }
                            Err(e) if e == "Cancelled" => {}
                            Err(e) => {
                                infer_failed.store(true, Ordering::SeqCst);
                                state.add_event(format!("  FAIL  {} — infer_experiment: {e}", sample.name));
                            }
                        }
                });

                // read_distribution.py
                let t_rdist = s.spawn(|| {
                    let script = config.rseqc_env.join("bin/read_distribution.py");
                        let script_str = match script.to_str() {
                            Some(s) => s.to_string(),
                            None => {
                                rdist_failed.store(true, Ordering::SeqCst);
                                state.add_event(format!("  FAIL  {} — read_distribution script path is not valid UTF-8", sample.name));
                                return;
                            }
                        };
                        let mut cmd = Command::new(&rseqc_python);
                        cmd.args([
                            script_str.as_str(),
                            "-i", bam_str.as_str(),
                            "-r", bed_str.as_str(),
                        ]);
                        match run_cancellable_capture(cmd) {
                            Ok(Some(stdout)) => {
                                if let Err(e) = atomic_write(&rdist_out, &stdout) {
                                    state.add_event(format!("  FAIL  {} — read_distribution write failed: {e}", sample.name));
                                    rdist_failed.store(true, Ordering::SeqCst);
                                } else {
                                    state.add_event(format!("  DONE  {} — read_distribution", sample.name));
                                }
                            }
                            Ok(None) => {
                                rdist_failed.store(true, Ordering::SeqCst);
                                state.add_event(format!("  FAIL  {} — read_distribution exited non-zero", sample.name));
                            }
                            Err(e) if e == "Cancelled" => {}
                            Err(e) => {
                                rdist_failed.store(true, Ordering::SeqCst);
                                state.add_event(format!("  FAIL  {} — read_distribution: {e}", sample.name));
                            }
                        }
                });

                // geneBody_coverage2.py (bigwig from Phase 2 is already present)
                let t_genebody = s.spawn(|| {
                    let script = config.rseqc_env.join("bin/geneBody_coverage2.py");
                        let script_str = match script.to_str() {
                            Some(s) => s.to_string(),
                            None => {
                                genebody_failed.store(true, Ordering::SeqCst);
                                state.add_event(format!("  FAIL  {} — geneBody_coverage2 script path is not valid UTF-8", sample.name));
                                return;
                            }
                        };
                        let log_dir = config.output_dir.join("logs");
                        let (gb_stdout, gb_stderr) = match make_log_stdio(&log_dir, &format!("{}.genebody", sample.name)) {
                            Ok(p) => p,
                            Err(e) => {
                                genebody_failed.store(true, Ordering::SeqCst);
                                state.add_event(format!("  FAIL  {} — geneBody_coverage2 log: {e}", sample.name));
                                return;
                            }
                        };
                        let mut cmd = Command::new(&rseqc_python);
                        cmd.args([
                            script_str.as_str(),
                            "-r", bed_str.as_str(),
                            "-i", bw_str.as_str(),
                            "-o", genebody_prefix.as_str(),
                        ]);
                        cmd.stdout(gb_stdout).stderr(gb_stderr);
                        match run_cancellable(cmd) {
                            Ok(true) => {
                                state.add_event(format!("  DONE  {} — geneBody_coverage2", sample.name));

                                // Generate ggplot2 PDF plots — required outputs
                                // included in the genebody SHA256 checkpoint.
                                {
                                    let genebody_prefix_str = genebody_prefix.to_string();
                                    let mut pdf_ok = true;

                                    // Curves plot
                                    match create_curves_r_script(&genebody_prefix_str) {
                                        Ok(_) => {
                                            let curves_r = format!("{}.geneBodyCoverage.curves.r", genebody_prefix_str);
                                            let curves_pdf = format!("{}.geneBodyCoverage.curves.pdf", genebody_prefix_str);
                                            match run_r_script(&curves_r) {
                                                Ok(_) => {
                                                    match validate_pdf_created(&curves_pdf) {
                                                        Ok(_) => state.add_event(format!("  DONE  {} — geneBodyCoverage curves plot", sample.name)),
                                                        Err(e) => {
                                                            state.add_event(format!("  FAIL  {} — curves PDF validation: {e}", sample.name));
                                                            pdf_ok = false;
                                                        }
                                                    }
                                                }
                                                Err(e) if e == "Cancelled" => {}
                                                Err(e) => {
                                                    state.add_event(format!("  FAIL  {} — geneBodyCoverage curves plot: {e}", sample.name));
                                                    pdf_ok = false;
                                                }
                                            }
                                        }
                                        Err(e) => {
                                            state.add_event(format!("  FAIL  {} — curves R script generation: {e}", sample.name));
                                            pdf_ok = false;
                                        }
                                    }

                                    // Heatmap plot (independent of curves)
                                    match create_heatmap_r_script(&genebody_prefix_str) {
                                        Ok(_) => {
                                            let heatmap_r = format!("{}.geneBodyCoverage.heatmap.r", genebody_prefix_str);
                                            let heatmap_pdf = format!("{}.geneBodyCoverage.heatMap.pdf", genebody_prefix_str);
                                            match run_r_script(&heatmap_r) {
                                                Ok(_) => {
                                                    match validate_pdf_created(&heatmap_pdf) {
                                                        Ok(_) => state.add_event(format!("  DONE  {} — geneBodyCoverage heatmap plot", sample.name)),
                                                        Err(e) => {
                                                            state.add_event(format!("  FAIL  {} — heatmap PDF validation: {e}", sample.name));
                                                            pdf_ok = false;
                                                        }
                                                    }
                                                }
                                                Err(e) if e == "Cancelled" => {}
                                                Err(e) => {
                                                    state.add_event(format!("  FAIL  {} — geneBodyCoverage heatmap plot: {e}", sample.name));
                                                    pdf_ok = false;
                                                }
                                            }
                                        }
                                        Err(e) => {
                                            state.add_event(format!("  FAIL  {} — heatmap R script generation: {e}", sample.name));
                                            pdf_ok = false;
                                        }
                                    }

                                    if !pdf_ok {
                                        genebody_failed.store(true, Ordering::SeqCst);
                                    }
                                }
                            }
                            Ok(false) => {
                                genebody_failed.store(true, Ordering::SeqCst);
                                state.add_event(format!("  FAIL  {} — geneBody_coverage2 exited non-zero", sample.name));
                            }
                            Err(e) if e == "Cancelled" => {}
                            Err(e) => {
                                genebody_failed.store(true, Ordering::SeqCst);
                                state.add_event(format!("  FAIL  {} — geneBody_coverage2: {e}", sample.name));
                            }
                        }
                });

                let _ = (t_infer.join(), t_rdist.join(), t_genebody.join());
            });
        }

        state.record_duration(t_start.elapsed().as_secs_f64());

        // On cancellation: clean up all partial outputs
        if is_cancelled() {
            let _ = fs::remove_file(&strand_out);
            let _ = fs::remove_file(&rdist_out);
            for suffix in &[
                ".geneBodyCoverage.txt", ".geneBodyCoverage_plot.r",
                ".geneBodyCoverage.pdf", ".geneBodyCoverage.curves.r",
                ".geneBodyCoverage.curves.pdf", ".geneBodyCoverage.heatmap.r",
                ".geneBodyCoverage.heatMap.pdf",
            ] {
                let _ = fs::remove_file(qc_dir.join(format!("{}{suffix}", sample.name)));
            }
            return Err("Cancelled".to_string());
        }

        // Per-tool checkpoint writes and cleanup.
        // Each tool that succeeded gets its checkpoint written immediately.
        // Each tool that failed gets only its own outputs cleaned — successful
        // sibling outputs are preserved for the next resume.
        let mut failed_tools: Vec<&str> = Vec::new();

        if infer_failed.load(Ordering::SeqCst) {
            let _ = fs::remove_file(&strand_out);
            remove_step_checkpoint(&config.output_dir, &sample.name, STEP_INFER);
            failed_tools.push("infer_experiment");
        } else {
            if let Err(e) = write_step_checkpoint(&config.output_dir, &sample.name, STEP_INFER) {
                state.add_event(format!("  FAIL  {} — infer checkpoint write: {e}", sample.name));
                failed_tools.push("infer_experiment(checkpoint)");
            }
        }

        if rdist_failed.load(Ordering::SeqCst) {
            let _ = fs::remove_file(&rdist_out);
            remove_step_checkpoint(&config.output_dir, &sample.name, STEP_RDIST);
            failed_tools.push("read_distribution");
        } else {
            if let Err(e) = write_step_checkpoint(&config.output_dir, &sample.name, STEP_RDIST) {
                state.add_event(format!("  FAIL  {} — rdist checkpoint write: {e}", sample.name));
                failed_tools.push("read_distribution(checkpoint)");
            }
        }

        if genebody_failed.load(Ordering::SeqCst) {
            for suffix in &[
                ".geneBodyCoverage.txt", ".geneBodyCoverage_plot.r",
                ".geneBodyCoverage.pdf", ".geneBodyCoverage.curves.r",
                ".geneBodyCoverage.curves.pdf", ".geneBodyCoverage.heatmap.r",
                ".geneBodyCoverage.heatMap.pdf",
            ] {
                let _ = fs::remove_file(qc_dir.join(format!("{}{suffix}", sample.name)));
            }
            remove_step_checkpoint(&config.output_dir, &sample.name, STEP_GENEBODY);
            failed_tools.push("geneBody_coverage2");
        } else {
            if let Err(e) = write_step_checkpoint(&config.output_dir, &sample.name, STEP_GENEBODY) {
                state.add_event(format!("  FAIL  {} — genebody checkpoint write: {e}", sample.name));
                failed_tools.push("geneBody_coverage2(checkpoint)");
            }
        }

        if !failed_tools.is_empty() {
            return Err(format!("{}: failed tools: {}", sample.name, failed_tools.join(", ")));
        }
    }

    Ok(())
}

// ─── Partial cleanup on failure ──────────────────────────────────────────────

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
            if let Err(e) = fs::remove_file(&path) {
                eprintln!("WARNING: could not remove partial STAR output {}: {e} \
                    — manual removal may be needed to avoid stale data", path.display());
            }
        }
    }
    // Remove _STARgenome, _STARpass1, and any _STARtmp* variant dirs
    // STAR may suffix _STARtmp with a PID or timestamp depending on version
    for subdir in ["_STARgenome", "_STARpass1"] {
        let dir = star_dir.join(format!("{sample_name}{subdir}"));
        if dir.exists() {
            let _ = fs::remove_dir_all(&dir);
        }
    }
    // Glob for any _STARtmp* directories (version-dependent suffixes).
    // SAFETY: sample_name is validated in discover_samples to [A-Za-z0-9_.-] so it
    // contains no glob metacharacters; the pattern is safe to pass to glob::glob.
    let tmp_pattern = star_dir.join(format!("{sample_name}_STARtmp*")).to_string_lossy().to_string();
    if let Ok(entries) = glob::glob(&tmp_pattern) {
        for entry in entries.flatten() {
            if entry.is_dir() {
                let _ = fs::remove_dir_all(&entry);
            }
        }
    }
}

/// Remove partial deeptools outputs for a sample (bigwig file).
fn cleanup_partial_deeptools(output_dir: &Path, sample_name: &str) {
    let bw_path = output_dir.join("bigwig").join(format!("{sample_name}.bw"));
    if bw_path.exists() {
        if let Err(e) = fs::remove_file(&bw_path) {
            eprintln!("WARNING: could not remove partial bigwig {}: {e}", bw_path.display());
        }
    }
}

/// Remove partial RSeQC/QC outputs for a sample (strand, read_dist, genebody + PDFs).
fn cleanup_partial_rseqc(output_dir: &Path, sample_name: &str) {
    let qc_dir = output_dir.join("qc");
    let suffixes = [
        ".strand.txt",
        ".read_distribution.txt",
        ".geneBodyCoverage.txt",
        ".geneBodyCoverage_plot.r",
        ".geneBodyCoverage.pdf",
        ".geneBodyCoverage.curves.r",
        ".geneBodyCoverage.curves.pdf",
        ".geneBodyCoverage.heatmap.r",
        ".geneBodyCoverage.heatMap.pdf",
    ];
    for suffix in &suffixes {
        let path = qc_dir.join(format!("{sample_name}{suffix}"));
        if path.exists() {
            if let Err(e) = fs::remove_file(&path) {
                eprintln!("WARNING: could not remove partial QC output {}: {e}", path.display());
            }
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
    if s.is_nan() || s.is_infinite() || s < 0.0 {
        return "??:??".to_string();
    }
    // Cap at 99:59:59 (359999 s) — beyond that, display "99:59:59+".
    // Avoids both the u64::MAX→f64 precision loss and the Duration nanosecond overflow
    // that occurs when secs_f64 × 1e9 exceeds u64::MAX.
    if s >= 359_999.0 {
        return "99:59:59+".to_string();
    }
    fmt_duration(Duration::from_secs_f64(s))
}

// ─── Layout helpers ──────────────────────────────────────────────────────────

fn truncate_to(s: &str, n: usize) -> String {
    if s.chars().count() > n {
        s.chars().take(n).collect()
    } else {
        s.to_string()
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum LayoutMode {
    Compact,
    Normal,
    Wide,
}

struct LayoutMetrics {
    overall_bar_w: usize,   // fills terminal: w - 2  (2 indent + 5 for " 100%")
    sample_bar_w: usize,    // fills terminal on Wide; capped on Compact/Normal
    max_name: usize,        // name column width, mode-dependent
    stats_rows: usize,      // 1 (Wide) or 2 (Compact/Normal)
}

fn compute_layout(w: usize, _h: usize) -> LayoutMetrics {
    let mode = match w {
        0..=79 => LayoutMode::Compact,
        80..=119 => LayoutMode::Normal,
        _ => LayoutMode::Wide,
    };

    // Overall bar: "  [bar] 100%"  →  2 + bar + 5 = w  →  bar = w - 7
    let overall_bar_w = w.saturating_sub(7).max(4);

    // Sample bar: "    [bar] 100%"  →  4 + bar + 5 = w  →  bar = w - 9
    // Compact: cap 28, Normal: cap 60, Wide: uncapped (fills terminal)
    let raw = w.saturating_sub(9);
    let sample_bar_w = match mode {
        LayoutMode::Compact => raw.min(28).max(4),
        LayoutMode::Normal => raw.min(60).max(4),
        LayoutMode::Wide => raw.max(4),
    };

    // Name column: leave budget for "  S " (4) + "[step_label] " (~14) + "HH:MM:SS / ~HH:MM:SS" (21) = 39
    let raw_name = w.saturating_sub(39);
    let max_name = match mode {
        LayoutMode::Compact => raw_name.min(16).max(6),
        LayoutMode::Normal => raw_name.min(28).max(8),
        LayoutMode::Wide => raw_name.min(48).max(10),
    };

    // Stats: 2 lines on Compact/Normal, 1 line on Wide
    let stats_rows = if mode == LayoutMode::Wide { 1 } else { 2 };

    LayoutMetrics {
        overall_bar_w,
        sample_bar_w,
        max_name,
        stats_rows,
    }
}

#[allow(unused_assignments)]
fn render_screen(
    stdout: &mut io::Stdout,
    state: &ProgressState,
    parallel_jobs: usize,
    resumed: usize,
) {
    let (term_w, term_h) = terminal::size().unwrap_or((80, 24));
    let w = term_w as usize;
    let h = term_h as usize;

    let lm = compute_layout(w, h);

    let elapsed = state.start_time.elapsed();
    let done = state.done_count();
    let total = state.total;
    let completed = state.completed.load(Ordering::Relaxed);
    let skipped = state.skipped.load(Ordering::Relaxed);
    let failed = state.failed.load(Ordering::Relaxed);
    let remaining = total.saturating_sub(done);
    let phase = state.phase();
    let avg_dur = state.avg_duration();

    let pct = if total > 0 { done.min(total) * 100 / total } else { 0 };
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

    let bar_width = lm.overall_bar_w;
    let filled = if total > 0 {
        bar_width * done.min(total) / total
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
        .map(|e| e.iter().cloned().collect())
        .unwrap_or_default();

    // Move cursor to top-left and clear from there
    let _ = execute!(
        stdout,
        cursor::MoveTo(0, 0),
        terminal::Clear(ClearType::All)
    );

    let mut row: u16 = 0;
    let footer_row = h.saturating_sub(1) as u16;  // Sticky footer at bottom

    // Top separator
    let _ = execute!(stdout, cursor::MoveTo(0, row), style::SetForegroundColor(Color::Cyan));
    println!("{}", "═".repeat(w));
    row += 1;

    // Header — Title (centered)
    let title = concat!("STAR-RSeQC v", env!("CARGO_PKG_VERSION"));
    let title_chars = title.chars().count();
    let title_pad_left = w.saturating_sub(title_chars) / 2;
    let title_pad_right = w.saturating_sub(title_chars).saturating_sub(title_pad_left);
    let title_line = format!("{}{}{}", " ".repeat(title_pad_left), title, " ".repeat(title_pad_right));
    let _ = execute!(stdout, cursor::MoveTo(0, row), style::SetForegroundColor(Color::White), style::SetAttribute(Attribute::Bold));
    println!("{}", title_line);
    row += 1;

    // Subtitle (centered, truncated if needed)
    let subtitle = "STAR 2-Pass Alignment + RSeQC Quality Control | Paired-End RNA-seq";
    let sub_display = truncate_to(subtitle, w);
    let sub_len = sub_display.chars().count();
    let sub_pad_left = w.saturating_sub(sub_len) / 2;
    let sub_pad_right = w.saturating_sub(sub_len).saturating_sub(sub_pad_left);
    let sub_line = format!("{}{}{}", " ".repeat(sub_pad_left), sub_display, " ".repeat(sub_pad_right));
    let _ = execute!(stdout, cursor::MoveTo(0, row), style::SetForegroundColor(Color::DarkGrey));
    println!("{}", sub_line);
    row += 1;

    // Separator
    let _ = execute!(stdout, cursor::MoveTo(0, row), style::SetForegroundColor(Color::Cyan));
    println!("{}", "═".repeat(w));
    row += 1;

    // Phase indicator
    let phase_line = format!("  {phase}");
    let _ = execute!(stdout, cursor::MoveTo(0, row), style::SetForegroundColor(Color::Magenta));
    println!("{}", truncate_to(&phase_line, w));
    row += 1;

    if resumed > 0 {
        let resume_line = format!(
            "  Resumed: {resumed} sample(s) already completed from previous run"
        );
        let _ = execute!(stdout, cursor::MoveTo(0, row), style::SetForegroundColor(Color::DarkYellow));
        println!("{}", truncate_to(&resume_line, w));
        row += 1;
    }

    // Separator
    let _ = execute!(stdout, cursor::MoveTo(0, row), style::SetForegroundColor(Color::Cyan));
    println!("{}", "═".repeat(w));
    row += 1;

    // Overall progress bar label
    let _ = execute!(stdout, cursor::MoveTo(0, row), style::SetForegroundColor(Color::White));
    println!("  OVERALL PROGRESS");
    row += 1;

    let _ = execute!(
        stdout,
        cursor::MoveTo(0, row),
        style::SetForegroundColor(Color::Green)
    );
    print!("  {bar_filled}");
    let _ = execute!(stdout, style::SetForegroundColor(Color::DarkGrey));
    print!("{bar_empty}");
    let _ = execute!(stdout, style::SetForegroundColor(Color::White));
    let pct_str = format!(" {:>3}%", pct);
    print!("{pct_str}");
    let _ = stdout.flush();
    row += 1;

    // Calculate estimated completion time
    let completion_time = if eta.as_secs() > 0 {
        let now = chrono::Local::now();
        let duration_secs = eta.as_secs() as i64;
        let completion = now + chrono::Duration::seconds(duration_secs);
        completion.format("%H:%M:%S").to_string()
    } else {
        "—".to_string()
    };

    if lm.stats_rows == 2 {
        // Line A: progress info
        let s = format!(
            "  {}/{} done   Elapsed: {}   ETA: {}",
            done, total, fmt_duration(elapsed), fmt_duration(eta)
        );
        let _ = execute!(stdout, cursor::MoveTo(0, row), style::SetForegroundColor(Color::White));
        println!("{}", truncate_to(&s, w));
        row += 1;
        // Line B: completion info
        let s = format!("  Complete: {}   Speed: {:.1}/min", completion_time, speed);
        let _ = execute!(stdout, cursor::MoveTo(0, row), style::SetForegroundColor(Color::White));
        println!("{}", truncate_to(&s, w));
        row += 1;
    } else {
        let s = format!(
            "  {}/{} done   Elapsed: {}   ETA: {}   Complete: {}   Speed: {:.1}/min",
            done, total, fmt_duration(elapsed), fmt_duration(eta), completion_time, speed
        );
        let _ = execute!(stdout, cursor::MoveTo(0, row), style::SetForegroundColor(Color::White));
        println!("{}", truncate_to(&s, w));
        row += 1;
    }

    let _ = execute!(stdout, cursor::MoveTo(0, row), style::SetForegroundColor(Color::Cyan));
    println!("{}", "═".repeat(w));
    row += 1;

    // Active jobs section
    let _ = execute!(stdout, cursor::MoveTo(0, row), style::SetForegroundColor(Color::Yellow));
    println!("  ACTIVE JOBS ({}/{})", active_count, parallel_jobs);
    row += 1;

    let active_jobs: Vec<(usize, &JobSlot)> = active_snapshot
        .iter()
        .enumerate()
        .filter_map(|(i, s)| s.as_ref().map(|j| (i, j)))
        .collect();

    // Guard against very small terminals (e.g. tmux panes): ensure row never exceeds
    // footer_row so MoveTo coordinates always stay within the terminal.
    let max_active_rows = if row as usize + 11 + lm.stats_rows < h {
        (h - row as usize - 11 - lm.stats_rows) / 2
    } else {
        0
    };
    let spinner_chars = ['|', '/', '-', '\\'];
    let spin_idx = (elapsed.as_millis() / 250) as usize;
    let sample_bar_w = lm.sample_bar_w;

    if active_jobs.is_empty() {
        let _ = execute!(stdout, cursor::MoveTo(0, row), style::SetForegroundColor(Color::DarkGrey));
        println!("  No active jobs");
        row += 1;
    }

    for (shown, (i, job)) in active_jobs.iter().enumerate() {
        if shown >= max_active_rows {
            let hidden = active_jobs.len().saturating_sub(shown);
            if hidden > 0 {
                let more = format!("  ... and {hidden} more active");
                let _more_len = more.chars().count();
                let _ = execute!(
                    stdout,
                    cursor::MoveTo(0, row),
                    style::SetForegroundColor(Color::DarkGrey)
                );
                print!("{}", more);
                row += 1;
            }
            break;
        }

        let spin = spinner_chars[(spin_idx + i) % 4];
        let job_elapsed_secs = job.started.elapsed().as_secs_f64();
        let job_elapsed_str = fmt_secs(job_elapsed_secs);
        let max_name = lm.max_name;
        let sample_len = job.sample.chars().count();
        let name = if sample_len > max_name {
            format!("{}...", job.sample.chars().take(max_name.saturating_sub(3)).collect::<String>())
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
        let line = format!("  {spin} {:<width$} {} {}", name, step_label, eta_part, width = max_name);
        let _ = execute!(stdout, cursor::MoveTo(0, row), style::SetForegroundColor(Color::White));
        println!("{}", truncate_to(&line, w));
        row += 1;

        // Row 2: per-sample progress bar
        let _ = execute!(stdout, cursor::MoveTo(0, row));
        if avg_dur > 0.0 {
            let frac = (job_elapsed_secs / avg_dur).min(1.0);
            let s_filled = (sample_bar_w as f64 * frac) as usize;
            let s_empty = sample_bar_w.saturating_sub(s_filled);
            let _ = execute!(stdout, style::SetForegroundColor(Color::Yellow));
            print!("    {}", "\u{2588}".repeat(s_filled));
            let _ = execute!(stdout, style::SetForegroundColor(Color::DarkGrey));
            print!("{}", "\u{2591}".repeat(s_empty));
            let _ = execute!(stdout, style::SetForegroundColor(Color::White));
            let s_pct = format!(" {:>3}%", (frac * 100.0) as usize);
            print!("{s_pct}");
        } else {
            // Indeterminate: pulse animation
            let pulse_pos = (spin_idx + i * 3) % (sample_bar_w + 4);
            print!("    ");
            for p in 0..sample_bar_w {
                if p >= pulse_pos.saturating_sub(2) && p <= pulse_pos {
                    let _ = execute!(stdout, style::SetForegroundColor(Color::Yellow));
                    print!("\u{2588}");
                } else {
                    let _ = execute!(stdout, style::SetForegroundColor(Color::DarkGrey));
                    print!("\u{2591}");
                }
            }
        }
        let _ = stdout.flush();
        row += 1;
    }

    let _ = execute!(stdout, cursor::MoveTo(0, row), style::SetForegroundColor(Color::Cyan));
    println!("{}", "═".repeat(w));
    row += 1;

    // Counters
    let _ = execute!(stdout, cursor::MoveTo(0, row));
    let _ = execute!(stdout, style::SetForegroundColor(Color::Green));
    print!("  Completed: {completed}");
    let _ = execute!(stdout, style::SetForegroundColor(Color::Yellow));
    print!("   Skipped: {skipped}");
    let _ = execute!(stdout, style::SetForegroundColor(Color::Red));
    print!("   Failed: {failed}");
    let _ = execute!(stdout, style::SetForegroundColor(Color::White));
    print!("   Remaining: {remaining}");
    let _ = stdout.flush();
    row += 1;

    let _ = execute!(stdout, cursor::MoveTo(0, row), style::SetForegroundColor(Color::Cyan));
    println!("{}", "═".repeat(w));
    row += 1;

    // Recent activity
    let _ = execute!(stdout, cursor::MoveTo(0, row), style::SetForegroundColor(Color::Magenta));
    println!("  RECENT ACTIVITY");
    row += 1;

    // Leave room for footer (2 rows: separator + footer)
    let max_event_rows = (footer_row as usize).saturating_sub(row as usize + 2).max(0);
    let start = events.len().saturating_sub(max_event_rows);
    for event_line in &events[start..] {
        let _ = execute!(stdout, cursor::MoveTo(0, row));
        let ev = truncate_to(event_line, w);
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
        let _ = stdout.flush();
        row += 1;
    }

    // === FOOTER (STICKY AT BOTTOM) ===
    // Footer separator (above footer)
    let _ = execute!(stdout, cursor::MoveTo(0, footer_row.saturating_sub(1)), style::SetForegroundColor(Color::Cyan));
    println!("{}", "═".repeat(w));

    // Footer info line (always at very bottom)
    let cancel_hint = if is_cancelled() {
        "CANCELLING..."
    } else {
        "Ctrl+C to cancel"
    };
    let timestamp = format!("Updated: {}", Local::now().format("%H:%M:%S"));
    let cancel_hint_len = cancel_hint.chars().count();
    let timestamp_len = timestamp.chars().count();
    let footer_pad = w
        .saturating_sub(cancel_hint_len)
        .saturating_sub(timestamp_len)
        .saturating_sub(2);

    let _ = execute!(stdout, cursor::MoveTo(0, footer_row));
    if is_cancelled() {
        let _ = execute!(stdout, style::SetForegroundColor(Color::Red));
    } else {
        let _ = execute!(stdout, style::SetForegroundColor(Color::DarkGrey));
    }
    print!("  {}", cancel_hint);
    let _ = execute!(stdout, style::SetForegroundColor(Color::DarkGrey));
    print!("{}", " ".repeat(footer_pad));
    print!("{}", timestamp);
    let _ = stdout.flush();

    let _ = execute!(stdout, style::ResetColor);
    let _ = stdout.flush();
}

// ─── Display thread ──────────────────────────────────────────────────────────

struct DisplayThread {
    flag: Arc<AtomicBool>,
    handle: Option<std::thread::JoinHandle<()>>,
    is_tty: bool,
}

impl DisplayThread {
    fn start(state: Arc<ProgressState>, parallel_jobs: usize, resumed: usize, is_tty: bool) -> Self {
        let flag = Arc::new(AtomicBool::new(false));
        let display_flag = Arc::clone(&flag);

        let handle = std::thread::spawn(move || {
            let mut out = io::stdout();
            loop {
                // Only write ANSI/TUI output when running on a real terminal.
                // On a pipe or log file we skip rendering to avoid garbage output,
                // but still poll for Ctrl+C events (harmless no-op on non-TTY).
                if is_tty {
                    render_screen(&mut out, &state, parallel_jobs, resumed);
                }

                if event::poll(REFRESH_INTERVAL).unwrap_or(false) {
                    if let Ok(event::Event::Key(key)) = event::read() {
                        if key.code == event::KeyCode::Char('c')
                            && key.modifiers.contains(event::KeyModifiers::CONTROL)
                        {
                            CANCELLED.store(true, Ordering::SeqCst);
                        }
                    }
                }

                if display_flag.load(Ordering::SeqCst) || is_cancelled() {
                    if is_tty {
                        render_screen(&mut out, &state, parallel_jobs, resumed);
                    }
                    break;
                }
            }
        });

        Self {
            flag,
            handle: Some(handle),
            is_tty,
        }
    }

    fn stop(&mut self) {
        self.flag.store(true, Ordering::SeqCst);
        if let Some(h) = self.handle.take() {
            let _ = h.join();
        }
        if self.is_tty {
            let mut out = io::stdout();
            let _ = terminal::disable_raw_mode();
            let _ = execute!(out, cursor::Show, terminal::LeaveAlternateScreen);
        }
    }
}

impl Drop for DisplayThread {
    fn drop(&mut self) {
        if self.handle.is_some() {
            self.stop();
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
            s.spawn(move || {
                loop {
                    if is_cancelled() {
                        break;
                    }
                    let idx = next.fetch_add(1, Ordering::Relaxed);
                    if idx >= items.len() {
                        break;
                    }
                    if let Err(e) = w(&items[idx], slot) {
                        if e != "Cancelled" {
                            errs.lock().unwrap_or_else(|e| e.into_inner()).push(e);
                        }
                    }
                    st.clear_slot(slot);
                    std::thread::yield_now();
                }
            });
        }
    });

    errors.into_inner().unwrap()
}

// ─── Environment validation ──────────────────────────────────────────────────

/// Resolve a binary that may be a bare name (looked up in $PATH) or an absolute path.
fn is_executable(p: &Path) -> bool {
    use std::os::unix::fs::PermissionsExt;
    p.metadata()
        .map(|m| m.is_file() && m.permissions().mode() & 0o111 != 0)
        .unwrap_or(false)
}

fn resolve_binary(bin: &Path) -> Option<PathBuf> {
    if bin.is_absolute() || bin.components().count() > 1 {
        // Explicit path — check directly
        if is_executable(bin) { Some(bin.to_path_buf()) } else { None }
    } else {
        // Bare name — search $PATH
        env::var_os("PATH")
            .map(|path_var| {
                env::split_paths(&path_var)
                    .map(|dir| dir.join(bin))
                    .find(|p| is_executable(p))
            })
            .flatten()
    }
}

/// Parse a version string like "1.12", "1.9", "1.21" into (major, minor).
fn parse_samtools_version(output: &str) -> Option<(u32, u32)> {
    // samtools --version first line: "samtools 1.12\n..." or "samtools 1.9\n..."
    let line = output.lines().next()?;
    let ver = line.strip_prefix("samtools ")?;
    let mut parts = ver.split('.');
    let major: u32 = parts.next()?.parse().ok()?;
    let minor: u32 = parts.next()?.trim_end_matches(|c: char| !c.is_ascii_digit()).parse().ok()?;
    Some((major, minor))
}

fn validate_environment(config: &mut Config) -> Result<(), String> {
    // samtools is always needed (BAM indexing runs regardless of --skip-alignment)
    if resolve_binary(&config.samtools).is_none() {
        return Err(format!(
            "samtools not found: {}\nInstall samtools or specify with --samtools <full-path>",
            config.samtools.display()
        ));
    }

    // The `-@` threads flag for `samtools index` requires samtools >= 1.12.
    // Earlier versions silently treat it as an unknown option and fail every BAM index.
    {
        let ver_out = Command::new(&config.samtools)
            .arg("--version")
            .output()
            .ok()
            .map(|o| String::from_utf8_lossy(&o.stdout).into_owned())
            .unwrap_or_default();
        match parse_samtools_version(&ver_out) {
            Some((maj, min)) if maj > 1 || (maj == 1 && min >= 12) => {}
            Some((maj, min)) => {
                return Err(format!(
                    "samtools {maj}.{min} detected — version >= 1.12 required for multi-threaded \
                     indexing (samtools index -@).\nUpgrade samtools or set --samtools to a newer binary."
                ));
            }
            None => {
                eprintln!("WARNING: Could not parse samtools version — assuming it supports -@ (>= 1.12)");
            }
        }
    }

    if !config.skip_alignment {
        let star_bin = config.star_env.join("bin/STAR");
        if !star_bin.exists() {
            return Err(format!(
                "STAR binary not found: {}\nIs the --star-env path correct?",
                star_bin.display()
            ));
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
    }

    {
        let rseqc_python = config.rseqc_env.join("bin/python");
        if !rseqc_python.exists() {
            return Err(format!(
                "RSeQC python not found: {}\nIs the --rseqc-env path correct?",
                rseqc_python.display()
            ));
        }

        for script in [
            "infer_experiment.py",
            "geneBody_coverage2.py",
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

        let bam_coverage = config.deeptools_env.join("bin/bamCoverage");
        if !bam_coverage.exists() {
            return Err(format!(
                "bamCoverage not found: {}\nIs the --deeptools-env path correct?",
                bam_coverage.display()
            ));
        }

        // Rscript is required for gene body coverage PDF generation (ggplot2).
        if resolve_binary(Path::new("Rscript")).is_none() {
            return Err(
                "Rscript not found in PATH — required for gene body coverage PDF plots.\n\
                 Install R: apt-get install r-base && Rscript -e 'install.packages(\"ggplot2\")'"
                    .to_string(),
            );
        }
    }

    // GTF is needed only when BED must be auto-generated
    if config.bed.is_none() {
        if !config.gtf.exists() {
            return Err(format!("GTF not found: {}", config.gtf.display()));
        }
    }

    if !config.fastq_dir.exists() {
        return Err(format!(
            "FASTQ directory not found: {}",
            config.fastq_dir.display()
        ));
    }

    // zcat is required by STAR's --readFilesCommand
    if !config.skip_alignment {
        if resolve_binary(Path::new("zcat")).is_none() {
            return Err("zcat not found in $PATH — required by STAR for gzip FASTQ decompression.\nInstall gzip or ensure zcat is in $PATH.".to_string());
        }
    }

    Ok(())
}

// ─── Audit trail ─────────────────────────────────────────────────────────────

#[derive(serde::Serialize, Clone)]
struct FastqAuditEntry {
    sample: String,
    r1: String,
    r1_sha256: String,
    r1_bytes: u64,
    r2: String,
    r2_sha256: String,
    r2_bytes: u64,
}

/// Hash all FASTQ inputs for audit. Designed to run in a background thread
/// concurrently with Phase 1 so startup is not blocked.
fn hash_fastq_inputs(samples: &[Sample]) -> Vec<FastqAuditEntry> {
    let streaming_hex = |p: &Path| -> String {
        sha256_file(p)
            .map(|b| b.iter().map(|x| format!("{x:02x}")).collect::<String>())
            .unwrap_or_else(|e| {
                eprintln!("WARNING: cannot hash {}: {e}", p.display());
                format!("UNREADABLE:{e}")
            })
    };
    samples.iter().map(|s| FastqAuditEntry {
        sample: s.name.clone(),
        r1: s.r1.display().to_string(),
        r1_sha256: streaming_hex(&s.r1),
        r1_bytes: fs::metadata(&s.r1).map(|m| m.len()).unwrap_or(0),
        r2: s.r2.display().to_string(),
        r2_sha256: streaming_hex(&s.r2),
        r2_bytes: fs::metadata(&s.r2).map(|m| m.len()).unwrap_or(0),
    }).collect()
}

/// Pre-computed reference file SHA256 hashes — computed in a background thread
/// started at pipeline launch so `write_run_info` does not block at run end.
struct RefHashes {
    genome_params_sha256: String,
    gtf_sha256: String,
    bed_sha256: String,
}

/// Write run_info.json — immutable audit record for this pipeline invocation.
/// `fastq_inputs` and `ref_hashes` are pre-computed by background threads to avoid
/// blocking pipeline startup or run-end summary on large file hashing.
fn write_run_info(
    output_dir: &Path,
    config: &Config,
    bed_path: &Path,
    run_timestamp: &str,
    sample_count: usize,
    fastq_inputs: Vec<FastqAuditEntry>,
    ref_hashes: RefHashes,
) {
    // Software versions
    let star_ver    = get_tool_version(&config.star_env.join("bin/STAR"), "--version");
    let samtools_ver = get_tool_version(&config.samtools, "--version");
    let bamcov_ver  = get_tool_version(&config.deeptools_env.join("bin/bamCoverage"), "--version");
    let python_ver  = get_tool_version(&config.rseqc_env.join("bin/python"), "--version");
    let rseqc_ver   = get_tool_version(&config.rseqc_env.join("bin/infer_experiment.py"), "-v");

    let RefHashes { genome_params_sha256, gtf_sha256, bed_sha256 } = ref_hashes;

    let command_line = std::env::args().collect::<Vec<_>>().join(" ");

    #[derive(serde::Serialize)]
    struct RefFile {
        path: String,
        sha256: String,
    }
    #[derive(serde::Serialize)]
    struct SoftwareVersions {
        star: String,
        samtools: String,
        bamcoverage: String,
        rseqc_python: String,
        rseqc_package: String,
    }
    #[derive(serde::Serialize)]
    struct References {
        genome_dir: String,
        genome_params_sha256: String,
        gtf: RefFile,
        bed12: RefFile,
    }
    #[derive(serde::Serialize)]
    struct RunInfo {
        pipeline_version: &'static str,
        run_timestamp: String,
        operator: String,
        operator_note: &'static str,
        command_line: String,
        sample_count: usize,
        software_versions: SoftwareVersions,
        references: References,
        inputs: Vec<FastqAuditEntry>,
    }

    let info = RunInfo {
        pipeline_version: env!("CARGO_PKG_VERSION"),
        run_timestamp: run_timestamp.to_string(),
        operator: config.operator.clone(),
        operator_note: if std::env::args().any(|a| a == "--operator") {
            "Identity provided via --operator flag"
        } else {
            "Identity sourced from $USER/$LOGNAME — not cryptographically verified"
        },
        command_line,
        sample_count,
        software_versions: SoftwareVersions {
            star: star_ver,
            samtools: samtools_ver,
            bamcoverage: bamcov_ver,
            rseqc_python: python_ver,
            rseqc_package: rseqc_ver,
        },
        references: References {
            genome_dir: config.genome_dir.display().to_string(),
            genome_params_sha256,
            gtf: RefFile {
                path: config.gtf.display().to_string(),
                sha256: gtf_sha256,
            },
            bed12: RefFile {
                path: bed_path.display().to_string(),
                sha256: bed_sha256,
            },
        },
        inputs: fastq_inputs,
    };

    match serde_json::to_string_pretty(&info) {
        Ok(json) => {
            // Write a timestamped file so every invocation (including resumes) is preserved.
            // The generic run_info_latest.json is a convenience copy of the most recent run.
            // Both are written atomically (tmp + rename) to prevent partial-write corruption.
            let ts_safe = run_timestamp.replace([':', '+', ' '], "-");
            // Include PID so concurrent or rapid-retry invocations produce distinct filenames.
            let pid = std::process::id();
            let ts_file = format!("run_info_{ts_safe}_pid{pid}.json");
            if let Err(e) = atomic_write(&output_dir.join(&ts_file), json.as_bytes()) {
                eprintln!("ERROR: Failed to write {ts_file}: {e}");
            }
            // Also write run_info_latest.json so downstream tooling has a stable name.
            // NOTE: This file is overwritten on each invocation (including resume). The
            // timestamped files above preserve the full history.
            if let Err(e) = atomic_write(&output_dir.join("run_info_latest.json"), json.as_bytes()) {
                eprintln!("ERROR: Failed to write run_info_latest.json: {e}");
            }
        }
        Err(e) => eprintln!("ERROR: Failed to serialise run_info: {e}"),
    }
}

// ─── Main ────────────────────────────────────────────────────────────────────

fn main() -> ExitCode {
    // Install a panic hook that restores terminal state before printing the panic message.
    // Without this, a panic while the TUI is active leaves the user's terminal in raw mode
    // with the alternate screen active (invisible cursor, no echo, garbled input).
    let default_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        // Best-effort terminal restoration — ignore errors.
        let _ = crossterm::terminal::disable_raw_mode();
        let _ = crossterm::execute!(
            io::stdout(),
            crossterm::cursor::Show,
            crossterm::terminal::LeaveAlternateScreen
        );
        default_hook(info);
    }));

    let mut config = match parse_args() {
        Ok(c) => c,
        Err(true)  => return ExitCode::SUCCESS,  // --help printed
        Err(false) => return ExitCode::FAILURE,  // argument error
    };

    let run_timestamp = Local::now().format("%Y-%m-%dT%H:%M:%S%z").to_string();
    eprintln!("star-rseqc v{} | Run with -h or --help for usage information", env!("CARGO_PKG_VERSION"));
    eprintln!("Run timestamp: {}", run_timestamp);
    eprintln!();

    // ── Validate environment ──
    if let Err(e) = validate_environment(&mut config) {
        eprintln!("Environment check failed:\n  {e}");
        return ExitCode::FAILURE;
    }
    eprintln!("Environment OK.");
    eprintln!(
        "Resources: {} job(s) x {} thread(s)/job, {:.1} GB BAM sort RAM{}",
        config.parallel_jobs, config.threads_per_sample,
        config.bam_sort_ram as f64 / 1e9,
        if config.resources_auto { " [auto-detected]" } else { " [manual]" }
    );

    // ── Discover samples ──
    let mut all_samples = discover_samples(&config.fastq_dir);
    if all_samples.is_empty() {
        eprintln!(
            "No paired-end samples found in {}\n\
             Expected files matching *_1P.fastq.gz with corresponding *_2P.fastq.gz",
            config.fastq_dir.display()
        );
        return ExitCode::FAILURE;
    }
    eprintln!("Discovered {} paired-end samples.", all_samples.len());

    // Sort by FASTQ R1 file size ascending (small samples first, for balanced batching)
    all_samples.sort_by_key(|s| fs::metadata(&s.r1).map(|m| m.len()).unwrap_or(0));

    // ── Create output structure ──
    for subdir in ["star", "qc", "logs", "bigwig"] {
        if let Err(e) = fs::create_dir_all(config.output_dir.join(subdir)) {
            eprintln!("Cannot create {}/{}: {}", config.output_dir.display(), subdir, e);
            return ExitCode::FAILURE;
        }
    }
    if let Err(e) = fs::create_dir_all(checkpoint_dir(&config.output_dir)) {
        eprintln!("Cannot create checkpoint directory: {e}");
        return ExitCode::FAILURE;
    }

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
            let gtf_bytes = fs::metadata(&config.gtf).map(|m| m.len()).unwrap_or(0);
            eprintln!("Converting GTF → BED12 ({:.1} GB)...", gtf_bytes as f64 / 1e9);
            if gtf_bytes > 1_500_000_000 {
                eprintln!("  NOTE: Large GTF — conversion will use 2–4 GB RAM. \
                    Use --bed to supply a pre-converted BED12 and skip this step.");
            }
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

    // ── Reference integrity check — warn if genome params are unreadable ──
    let genome_params = config.genome_dir.join("genomeParameters.txt");
    if !genome_params.exists() {
        eprintln!("WARNING: genomeParameters.txt not found in genome_dir — reference provenance cannot be verified");
    }

    // ── Resume detection (SHA256 verification with parallelization) ──
    eprintln!("Checking resume status (parallel SHA256 verification on {} samples)...", all_samples.len());
    let resume_start = Instant::now();
    let resume_results = check_resume_all_parallel(&config.output_dir, &all_samples);
    let resume_elapsed = resume_start.elapsed().as_secs_f64();
    let num_threads = resume_check_threads(all_samples.len());
    eprintln!("  ✓ SHA256 check completed in {:.2}s ({} threads, {:.1} samples/sec)",
        resume_elapsed,
        num_threads,
        all_samples.len() as f64 / resume_elapsed.max(0.01)
    );

    let mut already_done: usize = 0;
    let mut output_changed: usize = 0;
    let mut star_to_process: Vec<&Sample> = Vec::new();
    let mut deeptools_to_process: Vec<&Sample> = Vec::new();
    let mut rseqc_to_process: Vec<&Sample> = Vec::new();

    // Build a map of sample name → status from parallel results
    let status_map: HashMap<String, ResumeStatus> = resume_results
        .into_iter()
        .map(|(name, status)| (name, status))
        .collect();

    for s in &all_samples {
        let status = status_map.get(&s.name).cloned().unwrap_or(ResumeStatus::NotDone);
        match status {
            ResumeStatus::AllDone => {
                already_done += 1;
            }
            ResumeStatus::Phase1Changed => {
                if config.skip_alignment {
                    let bam = config.output_dir.join("star")
                        .join(format!("{}_Aligned.sortedByCoord.out.bam", s.name));
                    if bam.exists() {
                        eprintln!("  Phase 1 CHANGED: {} — --skip-alignment set; re-running deeptools + RSeQC on existing BAM", s.name);
                        cleanup_partial_deeptools(&config.output_dir, &s.name);
                        cleanup_partial_rseqc(&config.output_dir, &s.name);
                        deeptools_to_process.push(s);
                        rseqc_to_process.push(s);
                    } else {
                        eprintln!("  WARNING: {} — Phase 1 outputs changed and BAM is absent, but --skip-alignment set; sample will be skipped", s.name);
                    }
                } else {
                    eprintln!("  Phase 1 CHANGED: {} — cleaning partial outputs, will re-run STAR + deeptools + RSeQC", s.name);
                    cleanup_partial_star(&config.output_dir.join("star"), &s.name);
                    cleanup_partial_deeptools(&config.output_dir, &s.name);
                    cleanup_partial_rseqc(&config.output_dir, &s.name);
                    star_to_process.push(s);
                    deeptools_to_process.push(s);
                    rseqc_to_process.push(s);
                }
                output_changed += 1;
            }
            ResumeStatus::Phase2Changed => {
                eprintln!("  Phase 2 CHANGED: {} — cleaning partial outputs, will re-run deeptools + RSeQC", s.name);
                cleanup_partial_deeptools(&config.output_dir, &s.name);
                cleanup_partial_rseqc(&config.output_dir, &s.name);
                output_changed += 1;
                deeptools_to_process.push(s);
                rseqc_to_process.push(s);
            }
            ResumeStatus::Phase3Changed => {
                eprintln!("  Phase 3 CHANGED: {} — cleaning partial QC outputs, will re-run RSeQC only", s.name);
                cleanup_partial_rseqc(&config.output_dir, &s.name);
                output_changed += 1;
                rseqc_to_process.push(s);
            }
            ResumeStatus::NotDone => {
                // Clean any partial outputs from a previous incomplete run
                if !config.skip_alignment {
                    cleanup_partial_star(&config.output_dir.join("star"), &s.name);
                }
                cleanup_partial_deeptools(&config.output_dir, &s.name);
                cleanup_partial_rseqc(&config.output_dir, &s.name);

                if !config.skip_alignment {
                    star_to_process.push(s);
                }
                deeptools_to_process.push(s);
                rseqc_to_process.push(s);
            }
        }
    }

    if already_done > 0 {
        eprintln!(
            "Resuming: {already_done}/{} samples verified (SHA256 OK), {} STAR + {} deeptools + {} RSeQC to process.",
            all_samples.len(),
            star_to_process.len(),
            deeptools_to_process.len(),
            rseqc_to_process.len()
        );
    }
    if output_changed > 0 {
        eprintln!(
            "  {output_changed} sample(s) have corrupted/changed outputs — will re-process."
        );
    }

    // ── Dry run ── (before thread spawning — no I/O waste on dry-run)
    if config.dry_run {
        println!();
        println!("Dry run — {} samples discovered:\n", all_samples.len());
        println!("  {:<25} {:<50} {}", "SAMPLE", "R1", "STATUS");
        println!("  {}", "-".repeat(90));
        for s in &all_samples {
            let status = match status_map.get(&s.name).cloned().unwrap_or(ResumeStatus::NotDone) {
                ResumeStatus::AllDone => "✓ ALL COMPLETE",
                ResumeStatus::Phase1Changed => "Phase 1 CHANGED — re-run STAR+deeptools+RSeQC",
                ResumeStatus::Phase2Changed => "Phase 2 CHANGED — re-run deeptools+RSeQC",
                ResumeStatus::Phase3Changed => "Phase 3 CHANGED — re-run RSeQC only",
                ResumeStatus::NotDone => "PENDING",
            };
            println!("  {:<25} {:<50} {}", s.name, s.r1.display(), status);
        }
        println!();
        println!("Resource plan:");
        let src = if config.resources_auto { "auto" } else { "manual" };
        println!("  Threads per sample : {} ({})", config.threads_per_sample, src);
        println!("  BAM sort RAM       : {:.1} GB ({})", config.bam_sort_ram as f64 / 1e9, src);
        println!("  STAR jobs          : {} ({})", config.parallel_star_jobs,
            if config.parallel_star_jobs != config.parallel_jobs { "manual" } else { src });
        println!("  deeptools jobs     : {} ({})", config.parallel_deeptools_jobs,
            if config.parallel_deeptools_jobs != config.parallel_jobs { "manual" } else { src });
        println!("  RSeQC jobs         : {} ({})", config.parallel_rseqc_jobs,
            if config.parallel_rseqc_jobs != config.parallel_jobs { "manual" } else { src });
        println!("  Output: {}", config.output_dir.display());
        println!("  BED12:  {}", bed_path.display());
        return ExitCode::SUCCESS;  // No threads spawned; no audit written for dry-run
    }

    // ── Spawn background FASTQ hashing thread ──
    // Hashing all FASTQ inputs (potentially 100s of GB) runs concurrently with
    // Phase 1 so startup is not blocked. Joined at each exit point via finish_audit!.
    eprintln!("Spawning background FASTQ hash thread ({} files)...", all_samples.len() * 2);
    let samples_for_hash = all_samples.clone();
    let mut fastq_hash_thread: Option<std::thread::JoinHandle<Vec<FastqAuditEntry>>> =
        Some(std::thread::spawn(move || hash_fastq_inputs(&samples_for_hash)));

    // ── Spawn background reference-file hashing thread ──
    // Hashing the GTF (up to ~1 GB) in the background so write_run_info does not
    // block the final summary at run end.
    {
        let _genome_params_ref = &genome_params; // silence unused warning
    }
    let ref_genome_params = genome_params.clone();
    let ref_gtf = config.gtf.clone();
    let ref_bed = bed_path.clone();
    let mut ref_hashes_thread: Option<std::thread::JoinHandle<RefHashes>> =
        Some(std::thread::spawn(move || {
            let streaming_hash = |p: &Path| -> String {
                sha256_file(p)
                    .map(|b| b.iter().map(|x| format!("{x:02x}")).collect::<String>())
                    .unwrap_or_else(|e| {
                        eprintln!("WARNING: Cannot hash {} for audit trail: {e}", p.display());
                        format!("UNREADABLE:{e}")
                    })
            };
            RefHashes {
                genome_params_sha256: streaming_hash(&ref_genome_params),
                gtf_sha256: streaming_hash(&ref_gtf),
                bed_sha256: streaming_hash(&ref_bed),
            }
        }));

    // Helper: join background threads and write the audit record.
    // Uses Option so it can be called at multiple exit points; subsequent calls are no-ops.
    macro_rules! finish_audit {
        () => {{
            if let Some(handle) = fastq_hash_thread.take() {
                let fastq_inputs = handle.join().unwrap_or_default();
                let ref_hashes = ref_hashes_thread.take()
                    .and_then(|h| h.join().ok())
                    .unwrap_or_else(|| RefHashes {
                        genome_params_sha256: "UNCOMPUTED".to_string(),
                        gtf_sha256: "UNCOMPUTED".to_string(),
                        bed_sha256: "UNCOMPUTED".to_string(),
                    });
                write_run_info(
                    &config.output_dir, &config, &bed_path, &run_timestamp,
                    all_samples.len(), fastq_inputs, ref_hashes,
                );
                eprintln!("Run info written: {}/run_info_latest.json", config.output_dir.display());
            }
        }};
    }

    if star_to_process.is_empty() && deeptools_to_process.is_empty() && rseqc_to_process.is_empty() {
        eprintln!("All samples already completed. Refreshing summary files...");
        write_summary_files(&config.output_dir, &all_samples, false, &run_timestamp);
        finish_audit!();
        eprintln!("Done.");
        return ExitCode::SUCCESS;
    }

    // ── Phase-specific job counts ──
    let parallel_star_jobs     = config.parallel_star_jobs;
    let parallel_rseqc_jobs    = config.parallel_rseqc_jobs;
    let parallel_deeptools_jobs = config.parallel_deeptools_jobs;
    let pipeline_start = Instant::now();

    let mut stdout = io::stdout();
    let is_tty = stdout.is_tty();

    // ─────────────────────────────────────────────────────────────────────────
    // PHASE 1: STAR Alignment
    // ─────────────────────────────────────────────────────────────────────────

    let mut phase1_completed = 0usize;
    let mut phase1_failed = 0usize;

    let config_ref = &config;

    if !star_to_process.is_empty() {
        let phase1_label = format!(
            "Phase 1/3 — STAR alignment ({} samples)",
            star_to_process.len()
        );
        let state1 = Arc::new(ProgressState::new(
            star_to_process.len(),
            parallel_star_jobs,
            &phase1_label,
        ));
        state1.skipped.store(already_done, Ordering::Relaxed);

        if is_tty {
            let _ = execute!(stdout, terminal::EnterAlternateScreen, cursor::Hide);
            let _ = terminal::enable_raw_mode();
        }

        let mut display1 =
            DisplayThread::start(Arc::clone(&state1), parallel_star_jobs, already_done, is_tty);

        let bed_ref1 = &bed_path;
        let _star_errors = run_work_queue(&star_to_process, parallel_star_jobs, &state1, |sample, slot| {
            // Remove stale STAR checkpoint before running — forces recomputation
            remove_step_checkpoint(&config_ref.output_dir, &sample.name, STEP_STAR);
            // Also invalidate downstream checkpoints since STAR outputs will change
            remove_step_checkpoint(&config_ref.output_dir, &sample.name, STEP_DEEPTOOLS);
            remove_step_checkpoint(&config_ref.output_dir, &sample.name, STEP_INFER);
            remove_step_checkpoint(&config_ref.output_dir, &sample.name, STEP_RDIST);
            remove_step_checkpoint(&config_ref.output_dir, &sample.name, STEP_GENEBODY);

            let result = run_star_sample(sample, config_ref, bed_ref1, &state1, slot);
            match &result {
                Ok(()) => {
                    // Write per-step checkpoint immediately after STAR success
                    match write_step_checkpoint(&config_ref.output_dir, &sample.name, STEP_STAR) {
                        Ok(_) => {
                            state1.completed.fetch_add(1, Ordering::Relaxed);
                        }
                        Err(e) => {
                            state1.add_event(format!("  FAIL  {} — checkpoint write failed (disk full?): {e}", sample.name));
                            state1.failed.fetch_add(1, Ordering::Relaxed);
                        }
                    }
                }
                Err(e) if e == "Cancelled" => {
                    state1.add_event(format!("  STOP  {} — cancelled", sample.name));
                }
                Err(e) => {
                    // STAR failed — no checkpoint written (step stays incomplete)
                    state1.failed.fetch_add(1, Ordering::Relaxed);
                    state1.add_event(format!("  FAIL  {} — {}", sample.name, e));
                }
            }
            result
        });

        display1.stop();

        phase1_completed = state1.completed.load(Ordering::Relaxed);
        phase1_failed = state1.failed.load(Ordering::Relaxed);
        let phase1_elapsed = fmt_duration(state1.start_time.elapsed());

        println!();
        println!("  ┌─────────────────────────────────────────────────────────┐");
        println!("  │             Phase 1 — STAR Alignment Complete            │");
        println!("  ├─────────────────────────────────────────────────────────┤");
        println!("  │  Completed:  {:<40}│", phase1_completed);
        println!("  │  Failed:     {:<40}│", phase1_failed);
        println!("  │  Elapsed:    {:<40}│", phase1_elapsed);
        println!("  └─────────────────────────────────────────────────────────┘");
        println!();

        if is_cancelled() {
            write_summary_files(&config.output_dir, &all_samples, true, &run_timestamp);
            println!();
            println!("  Pipeline cancelled by user during Phase 1.");
            println!();
            finish_audit!();
            return ExitCode::FAILURE;
        }
    }

    // ─────────────────────────────────────────────────────────────────────────
    // PHASE 2: deeptools — BAM → bigwig (bamCoverage)
    // ─────────────────────────────────────────────────────────────────────────

    // Only include samples where STAR BAM exists
    let deeptools_to_process_phase2: Vec<&Sample> = deeptools_to_process.iter()
        .filter(|s| {
            config.output_dir.join("star")
                .join(format!("{}_Aligned.sortedByCoord.out.bam", s.name))
                .exists()
        })
        .copied()
        .collect();

    let mut phase2_completed = 0;
    let mut phase2_failed_final = 0;

    if deeptools_to_process_phase2.is_empty() {
        println!("  No samples with completed STAR alignment for deeptools processing.");
    } else {
        let phase2_label = format!(
            "Phase 2/3 — deeptools bamCoverage ({} samples)",
            deeptools_to_process_phase2.len()
        );
        let state2 = Arc::new(ProgressState::new(
            deeptools_to_process_phase2.len(),
            parallel_deeptools_jobs,
            &phase2_label,
        ));

        if is_tty {
            let _ = execute!(stdout, terminal::EnterAlternateScreen, cursor::Hide);
            let _ = terminal::enable_raw_mode();
        }

        let mut display2 = DisplayThread::start(Arc::clone(&state2), parallel_deeptools_jobs, 0, is_tty);

        let _ = run_work_queue(&deeptools_to_process_phase2, parallel_deeptools_jobs, &state2, |sample, slot| {
            // Remove stale deeptools checkpoint before running
            remove_step_checkpoint(&config_ref.output_dir, &sample.name, STEP_DEEPTOOLS);
            // Invalidate downstream RSeQC checkpoints since bigwig will change
            remove_step_checkpoint(&config_ref.output_dir, &sample.name, STEP_INFER);
            remove_step_checkpoint(&config_ref.output_dir, &sample.name, STEP_RDIST);
            remove_step_checkpoint(&config_ref.output_dir, &sample.name, STEP_GENEBODY);

            let result = run_deeptools_phase(sample, config_ref, &state2, slot);
            match &result {
                Ok(()) => {
                    // Write per-step checkpoint immediately after deeptools success
                    match write_step_checkpoint(&config_ref.output_dir, &sample.name, STEP_DEEPTOOLS) {
                        Ok(_) => {
                            state2.completed.fetch_add(1, Ordering::Relaxed);
                        }
                        Err(e) => {
                            state2.add_event(format!("  FAIL  {} — checkpoint write failed (disk full?): {e}", sample.name));
                            state2.failed.fetch_add(1, Ordering::Relaxed);
                        }
                    }
                }
                Err(e) if e != "Cancelled" => {
                    // deeptools failed — no checkpoint written; STAR checkpoint is preserved
                    state2.failed.fetch_add(1, Ordering::Relaxed);
                }
                _ => {}
            }
            result
        });

        display2.stop();

        let phase2_total_elapsed = fmt_duration(state2.start_time.elapsed());
        phase2_completed = state2.completed.load(Ordering::Relaxed);
        phase2_failed_final = state2.failed.load(Ordering::Relaxed);

        println!();
        println!("  ┌──────────────────────────────────────────────────────────┐");
        println!("  │  Phase 2 — deeptools bamCoverage (BAM → bigwig)          │");
        println!("  ├──────────────────────────────────────────────────────────┤");
        println!("  │  Completed:  {:<41}│", phase2_completed);
        println!("  │  Failed:     {:<41}│", phase2_failed_final);
        println!("  │  Elapsed:    {:<41}│", phase2_total_elapsed);
        println!("  └──────────────────────────────────────────────────────────┘");
        println!();
    }

    if is_cancelled() {
        write_summary_files(&config.output_dir, &all_samples, true, &run_timestamp);
        println!();
        println!("  Pipeline cancelled by user during Phase 2.");
        println!();
        finish_audit!();
        return ExitCode::FAILURE;
    }

    // ─────────────────────────────────────────────────────────────────────────
    // PHASE 3: RSeQC — infer_experiment + read_distribution + geneBody_coverage2
    // ─────────────────────────────────────────────────────────────────────────

    // Only include samples where bigwig exists (Phase 2 complete)
    let rseqc_to_process_phase3: Vec<&Sample> = rseqc_to_process.iter()
        .filter(|s| {
            config.output_dir.join("bigwig")
                .join(format!("{}.bw", s.name))
                .exists()
        })
        .copied()
        .collect();

    let mut phase3_completed = 0;
    let mut phase3_failed_final = 0;

    if rseqc_to_process_phase3.is_empty() {
        println!("  No samples with completed bigwig for RSeQC processing.");
    } else {
        let bed_ref = &bed_path;

        let phase3_label = format!(
            "Phase 3/3 — RSeQC (infer + read_dist + geneBody) ({} samples)",
            rseqc_to_process_phase3.len()
        );
        let state3 = Arc::new(ProgressState::new(
            rseqc_to_process_phase3.len(),
            parallel_rseqc_jobs,
            &phase3_label,
        ));

        if is_tty {
            let _ = execute!(stdout, terminal::EnterAlternateScreen, cursor::Hide);
            let _ = terminal::enable_raw_mode();
        }

        let mut display3 = DisplayThread::start(Arc::clone(&state3), parallel_rseqc_jobs, 0, is_tty);

        let _ = run_work_queue(&rseqc_to_process_phase3, parallel_rseqc_jobs, &state3, |sample, slot| {
            let result = run_rseqc_phase3(sample, config_ref, bed_ref, &state3, slot);
            match &result {
                Ok(()) => {
                    // Checkpoints already written inside run_rseqc_phase3 per sub-tool
                    state3.completed.fetch_add(1, Ordering::Relaxed);
                }
                Err(e) if e != "Cancelled" => {
                    // Partial checkpoints written inside run_rseqc_phase3 for tools
                    // that succeeded; failed tools had their outputs + checkpoints cleaned.
                    state3.failed.fetch_add(1, Ordering::Relaxed);
                }
                _ => {}
            }
            result
        });

        display3.stop();

        let phase3_total_elapsed = fmt_duration(state3.start_time.elapsed());
        phase3_completed = state3.completed.load(Ordering::Relaxed);
        phase3_failed_final = state3.failed.load(Ordering::Relaxed);

        println!();
        println!("  ┌──────────────────────────────────────────────────────────┐");
        println!("  │  Phase 3 — RSeQC (infer + read_dist + geneBody_cov2)    │");
        println!("  ├──────────────────────────────────────────────────────────┤");
        println!("  │  Completed:  {:<41}│", phase3_completed);
        println!("  │  Failed:     {:<41}│", phase3_failed_final);
        println!("  │  Elapsed:    {:<41}│", phase3_total_elapsed);
        println!("  └──────────────────────────────────────────────────────────┘");
        println!();
    }

    let was_cancelled = is_cancelled();
    let total = all_samples.len();
    let elapsed_str = fmt_duration(pipeline_start.elapsed());

    // ── Write audit trail (join background hash thread) + summary files ──
    finish_audit!();
    write_summary_files(&config.output_dir, &all_samples, was_cancelled, &run_timestamp);

    // ── Final summary ──
    println!();
    println!("  \u{2554}{}\u{2557}", "\u{2550}".repeat(52));
    if was_cancelled {
        println!(
            "  \u{2551}         STAR-RSeQC  -  Cancelled by user           \u{2551}"
        );
    } else {
        println!(
            "  \u{2551}           STAR-RSeQC  -  Run Complete              \u{2551}"
        );
    }
    println!("  \u{2560}{}\u{2563}", "\u{2550}".repeat(52));
    println!("  \u{2551}  Total samples:      {:<29}\u{2551}", total);
    println!(
        "  \u{2551}  Phase 1 (STAR) complete:     {:<19}\u{2551}",
        phase1_completed
    );
    println!(
        "  \u{2551}  Phase 1 (STAR) failed:       {:<19}\u{2551}",
        phase1_failed
    );
    println!(
        "  \u{2551}  Phase 2 (deeptools) complete:{:<19}\u{2551}",
        phase2_completed
    );
    println!(
        "  \u{2551}  Phase 2 (deeptools) failed:  {:<19}\u{2551}",
        phase2_failed_final
    );
    println!(
        "  \u{2551}  Phase 3 (RSeQC) complete:    {:<19}\u{2551}",
        phase3_completed
    );
    println!(
        "  \u{2551}  Phase 3 (RSeQC) failed:      {:<19}\u{2551}",
        phase3_failed_final
    );
    if already_done > 0 {
        println!(
            "  \u{2551}  Resumed (SHA256 OK): {:<28}\u{2551}",
            already_done
        );
    }
    println!("  \u{2551}  Total elapsed:      {:<29}\u{2551}", elapsed_str);
    println!("  \u{2551}  Threads/sample:     {:<29}\u{2551}",
        if config.resources_auto { format!("{} (auto)", config.threads_per_sample) }
        else { config.threads_per_sample.to_string() });
    println!("  \u{2551}  STAR jobs:          {:<29}\u{2551}",
        if config.resources_auto { format!("{} (auto)", parallel_star_jobs) }
        else { parallel_star_jobs.to_string() });
    println!("  \u{2551}  deeptools jobs:     {:<29}\u{2551}",
        if config.resources_auto { format!("{} (auto)", parallel_deeptools_jobs) }
        else { parallel_deeptools_jobs.to_string() });
    println!("  \u{2551}  RSeQC jobs:         {:<29}\u{2551}",
        if config.resources_auto { format!("{} (auto)", parallel_rseqc_jobs) }
        else { parallel_rseqc_jobs.to_string() });
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

    let total_failures = phase1_failed + phase2_failed_final + phase3_failed_final;
    if total_failures > 0 {
        println!("  {} sample(s) failed during processing.", total_failures);
        println!();
    }

    if was_cancelled {
        println!("  Run was cancelled. Re-run the same command to resume.");
        println!();
    } else if total_failures > 0 {
        println!("  Some samples failed. Re-run to retry failed samples.");
        println!();
    } else if config.clean_sam {
        let cleaned = cleanup_sam_files(&config.output_dir, &all_samples);
        println!("  --clean-sam: deleted {} chimeric SAM file(s).", cleaned);
        println!();
    }

    if total_failures == 0 && !was_cancelled {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}

// ─── SAM file cleanup ────────────────────────────────────────────────────────

fn cleanup_sam_files(output_dir: &Path, samples: &[Sample]) -> usize {
    let star_dir = output_dir.join("star");
    let mut cleaned = 0;

    for sample in samples {
        let sam_path = star_dir.join(format!("{}_Chimeric.out.sam", sample.name));
        if sam_path.exists() {
            if let Ok(_) = fs::remove_file(&sam_path) {
                cleaned += 1;
            }
        }
    }

    cleaned
}

// ─── Summary files ───────────────────────────────────────────────────────────

fn write_summary_files(output_dir: &Path, samples: &[Sample], was_cancelled: bool, run_timestamp: &str) {
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

    let rows: Vec<SummaryRow> = samples
        .iter()
        .map(|s| {
            let n = &s.name;
            let digests = SampleDigests {
                star:      read_step_checkpoint(&output_dir, n, STEP_STAR).unwrap_or_else(|| "INCOMPLETE".to_string()),
                deeptools: read_step_checkpoint(&output_dir, n, STEP_DEEPTOOLS).unwrap_or_else(|| "INCOMPLETE".to_string()),
                infer:     read_step_checkpoint(&output_dir, n, STEP_INFER).unwrap_or_else(|| "INCOMPLETE".to_string()),
                rdist:     read_step_checkpoint(&output_dir, n, STEP_RDIST).unwrap_or_else(|| "INCOMPLETE".to_string()),
                genebody:  read_step_checkpoint(&output_dir, n, STEP_GENEBODY).unwrap_or_else(|| "INCOMPLETE".to_string()),
            };
            let sha256 = format!(
                "star:{}|deeptools:{}|infer:{}|rdist:{}|genebody:{}",
                digests.star,
                digests.deeptools,
                digests.infer,
                digests.rdist,
                digests.genebody,
            );
            SummaryRow {
                sample: n.clone(),
                sha256,
                // STAR
                log_final: star_dir.join(format!("{n}_Log.final.out")).exists(),
                log_out: star_dir.join(format!("{n}_Log.out")).exists(),
                log_progress: star_dir.join(format!("{n}_Log.progress.out")).exists(),
                bam_sorted: star_dir.join(format!("{n}_Aligned.sortedByCoord.out.bam")).exists(),
                bam_index: star_dir.join(format!("{n}_Aligned.sortedByCoord.out.bam.bai")).exists(),
                bam_transcriptome: star_dir.join(format!("{n}_Aligned.toTranscriptome.out.bam")).exists(),
                gene_counts: star_dir.join(format!("{n}_ReadsPerGene.out.tab")).exists(),
                splice_junctions: star_dir.join(format!("{n}_SJ.out.tab")).exists(),
                chimeric_junction: star_dir.join(format!("{n}_Chimeric.out.junction")).exists(),
                chimeric_sam: star_dir.join(format!("{n}_Chimeric.out.sam")).exists(),
                // RSeQC
                strand_qc: qc_dir.join(format!("{n}.strand.txt")).exists(),
                genebody_txt: qc_dir.join(format!("{n}.geneBodyCoverage.txt")).exists(),
                genebody_r: qc_dir.join(format!("{n}.geneBodyCoverage_plot.r")).exists(),
                genebody_curves_pdf: qc_dir.join(format!("{n}.geneBodyCoverage.curves.pdf")).exists(),
                genebody_heatmap_pdf: qc_dir.join(format!("{n}.geneBodyCoverage.heatMap.pdf")).exists(),
                readdist_qc: qc_dir.join(format!("{n}.read_distribution.txt")).exists(),
            }
        })
        .collect();

    // JSON — wrap rows in a root object with run_cancelled field
    #[derive(serde::Serialize)]
    struct Summary<'a> {
        run_cancelled: bool,
        samples: &'a Vec<SummaryRow>,
    }
    let summary = Summary { run_cancelled: was_cancelled, samples: &rows };
    match serde_json::to_string_pretty(&summary) {
        Ok(json) => {
            if let Err(e) = atomic_write(&output_dir.join("pipeline_summary.json"), json.as_bytes()) {
                eprintln!("ERROR: Failed to write pipeline_summary.json: {e}");
            }
        }
        Err(e) => eprintln!("ERROR: Failed to serialise pipeline_summary.json: {e}"),
    }

    // TSV (compact: group STAR and QC status)
    let mut tsv = String::new();
    tsv.push_str(&format!("# run_timestamp: {}\n", run_timestamp));
    if was_cancelled {
        tsv.push_str("# WARNING: Run was cancelled — outputs may be incomplete\n");
    }
    tsv.push_str(
        "sample\tsha256\trun_cancelled\tlog_final\tlog_out\tlog_progress\tbam_sorted\tbam_index\t\
         bam_transcriptome\tgene_counts\tsplice_junctions\tchimeric_junction\tchimeric_sam\t\
         strand_qc\tgenebody_txt\tgenebody_r\tgenebody_curves_pdf\tgenebody_heatmap_pdf\treaddist_qc\n",
    );
    for r in &rows {
        let ok = |b: bool| if b { "OK" } else { "MISSING" };
        tsv.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            r.sample, r.sha256,
            was_cancelled,
            ok(r.log_final), ok(r.log_out), ok(r.log_progress),
            ok(r.bam_sorted), ok(r.bam_index),
            ok(r.bam_transcriptome), ok(r.gene_counts), ok(r.splice_junctions),
            ok(r.chimeric_junction), ok(r.chimeric_sam), ok(r.strand_qc),
            ok(r.genebody_txt), ok(r.genebody_r), ok(r.genebody_curves_pdf),
            ok(r.genebody_heatmap_pdf), ok(r.readdist_qc),
        ));
    }
    if let Err(e) = atomic_write(&output_dir.join("pipeline_summary.tsv"), tsv.as_bytes()) {
        eprintln!("ERROR: Failed to write pipeline_summary.tsv: {e}");
    }
}
