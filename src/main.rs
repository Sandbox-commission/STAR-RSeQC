mod checkpoint;
mod config;
mod gtf;
mod pipeline;
mod sample;
mod tui;

use crossterm::{cursor, execute, terminal};
use log::{error, info, warn};
use std::fs;
use std::io;
use std::process::ExitCode;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use std::time::Duration;

use checkpoint::{check_resume, checkpoint_dir, ResumeStatus};
use config::VERSION;
use pipeline::{process_sample, run_work_queue, validate_environment, write_summary_files};
use sample::discover_samples;
use tui::{DisplayThread, ProgressState};

pub(crate) static CANCELLED: AtomicBool = AtomicBool::new(false);

pub(crate) fn is_cancelled() -> bool {
    CANCELLED.load(Ordering::Relaxed)
}

fn main() -> ExitCode {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info"))
        .format_timestamp(None)
        .format_target(false)
        .init();

    let config = match config::parse_args() {
        Some(c) => c,
        None => return ExitCode::SUCCESS,
    };

    info!(
        "star-rseqc v{} | Run with -h or --help for usage information",
        VERSION
    );

    // ── Validate environment ──
    if let Err(e) = validate_environment(&config) {
        error!("Environment check failed:\n  {e}");
        return ExitCode::FAILURE;
    }
    info!("Environment OK.");
    info!(
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
    let all_samples = discover_samples(&config.fastq_dir, &config.r1_suffix, &config.r2_suffix);
    if all_samples.is_empty() {
        error!(
            "No paired-end samples found in {}\n\
             Expected files matching *{}.fastq.gz with corresponding *{}.fastq.gz",
            config.fastq_dir.display(),
            config.r1_suffix,
            config.r2_suffix,
        );
        return ExitCode::FAILURE;
    }
    info!("Discovered {} paired-end samples.", all_samples.len());

    // ── Create output structure ──
    for subdir in ["star", "qc", "logs"] {
        if let Err(e) = fs::create_dir_all(config.output_dir.join(subdir)) {
            error!(
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
            error!("BED file not found: {}", bed.display());
            return ExitCode::FAILURE;
        }
        bed.clone()
    } else {
        let auto_bed = config.output_dir.join("annotation.bed12");
        if auto_bed.exists() {
            info!("Reusing cached BED12: {}", auto_bed.display());
        } else {
            info!("Converting GTF -> BED12...");
            match gtf::gtf_to_bed12(&config.gtf, &auto_bed) {
                Ok(n) => info!("BED12: {} transcripts written.", n),
                Err(e) => {
                    error!("GTF->BED12 failed: {e}");
                    return ExitCode::FAILURE;
                }
            }
        }
        auto_bed
    };

    // ── Resume detection (SHA256 verification) ──
    info!("Checking resume status (SHA256 verification)...");
    let mut already_done: usize = 0;
    let mut output_changed: usize = 0;
    let mut to_process: Vec<&sample::Sample> = Vec::new();

    for s in &all_samples {
        match check_resume(&config.output_dir, &s.name) {
            ResumeStatus::SameHash => {
                already_done += 1;
            }
            ResumeStatus::StarChanged { old, new } => {
                warn!(
                    "STAR changed: {} (was {}..., now {}...)",
                    s.name,
                    &old[..12.min(old.len())],
                    &new[..12.min(new.len())]
                );
                output_changed += 1;
                to_process.push(s);
            }
            ResumeStatus::RseqcChanged { old, new } => {
                warn!(
                    "RSeQC changed: {} (was {}..., now {}...)",
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
                warn!(
                    "STAR+RSeQC changed: {} (star: {}...->{}..., rseqc: {}...->{}...)",
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
        info!(
            "Resuming: {already_done}/{} samples verified (SHA256 OK), {} to process.",
            all_samples.len(),
            to_process.len()
        );
    }
    if output_changed > 0 {
        warn!("{output_changed} sample(s) have corrupted/changed outputs — will re-process.");
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
        if config.r1_suffix != "_1P" || config.r2_suffix != "_2P" {
            println!(
                "  FASTQ suffixes: R1='{}', R2='{}'",
                config.r1_suffix, config.r2_suffix
            );
        }
        if !config.star_extra_args.is_empty() {
            println!("  Extra STAR args: {}", config.star_extra_args.join(" "));
        }
        return ExitCode::SUCCESS;
    }

    if to_process.is_empty() {
        info!("All samples already completed. Nothing to do.");
        return ExitCode::SUCCESS;
    }

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
    state
        .skipped
        .store(already_done, Ordering::Relaxed);
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
                let digests =
                    checkpoint::sha256_outputs(&config_ref.output_dir, &sample.name);
                checkpoint::write_checkpoint(&config_ref.output_dir, &sample.name, &digests);
                state
                    .completed
                    .fetch_add(1, Ordering::Relaxed);
                state.add_event(format!(
                    "  DONE  {} — star:{}... rseqc:{}...",
                    sample.name,
                    &digests.star[..12],
                    &digests.rseqc[..12]
                ));
            }
            Err(e) if e == "Cancelled" => {
                state.add_event(format!("  STOP  {} — cancelled", sample.name));
            }
            Err(e) => {
                checkpoint::remove_checkpoint(&config_ref.output_dir, &sample.name);
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
    let elapsed_str = tui::fmt_duration_pub(state.start_time.elapsed());

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
