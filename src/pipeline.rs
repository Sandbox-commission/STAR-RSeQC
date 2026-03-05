use std::fs::{self, File};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use log::{info, warn};

use crate::checkpoint::{checkpoint_dir, parse_checkpoint};
use crate::config::Config;
use crate::sample::Sample;
use crate::tui::ProgressState;

// ─── Run command with cancellation ───────────────────────────────────────────

pub(crate) fn run_cancellable(mut cmd: Command, program_name: &str) -> Result<bool, String> {
    let mut child = cmd.spawn().map_err(|e| {
        format!(
            "Failed to launch '{program_name}': {e}. \
             Check that the binary exists and is executable."
        )
    })?;

    loop {
        match child.try_wait() {
            Ok(Some(status)) => return Ok(status.success()),
            Ok(None) => {
                if crate::is_cancelled() {
                    let _ = child.kill();
                    let _ = child.wait();
                    return Err("Cancelled".to_string());
                }
                std::thread::sleep(Duration::from_millis(200));
            }
            Err(e) => return Err(format!("Wait error for '{program_name}': {e}")),
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

pub(crate) fn process_sample(
    sample: &Sample,
    config: &Config,
    bed_path: &Path,
    state: &ProgressState,
    slot: usize,
) -> Result<(), String> {
    if crate::is_cancelled() {
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
            state.add_event(format!(
                "  RESUME  {} — BAM exists, skipping STAR",
                sample.name
            ));
        } else {
            let star_bin = config.star_env.join("bin/STAR");
            let log_name = format!("{}.star", sample.name);
            let (stdout_cfg, stderr_cfg) = make_log_stdio(&log_dir, &log_name);

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
            ]);

            // Append user-supplied extra STAR args
            if !config.star_extra_args.is_empty() {
                cmd.args(&config.star_extra_args);
            }

            cmd.stdout(stdout_cfg).stderr(stderr_cfg);

            let log_path = log_dir.join(format!("{log_name}.log"));
            match run_cancellable(cmd, "STAR") {
                Ok(true) => {
                    state.add_event(format!("  DONE  {} — STAR alignment", sample.name));
                }
                Ok(false) => {
                    cleanup_partial_star(&star_dir, &sample.name);
                    state.add_event(format!(
                        "  FAIL  {} — STAR alignment error (cleaned)",
                        sample.name
                    ));
                    return Err(format!(
                        "{}: STAR failed — check log at {}",
                        sample.name,
                        log_path.display()
                    ));
                }
                Err(e) => {
                    cleanup_partial_star(&star_dir, &sample.name);
                    return Err(e);
                }
            }
        }
    }

    if crate::is_cancelled() {
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

        match run_cancellable(cmd, "samtools") {
            Ok(true) => {
                state.add_event(format!("  DONE  {} — samtools index", sample.name));
            }
            Ok(false) => {
                state.add_event(format!("  FAIL  {} — samtools index error", sample.name));
                return Err(format!(
                    "{}: samtools index failed — verify samtools version (need 1.15+): {}",
                    sample.name,
                    config.samtools.display()
                ));
            }
            Err(e) => return Err(e),
        }
    }

    if crate::is_cancelled() {
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
            let output = cmd.output().map_err(|e| {
                format!("{}: Failed to run {}: {}", sample.name, script.display(), e)
            })?;
            if output.status.success() {
                fs::write(&strand_out, &output.stdout)
                    .map_err(|e| format!("Write {}: {}", strand_out.display(), e))?;
                state.add_event(format!("  DONE  {} — infer_experiment", sample.name));
            } else {
                let stderr = String::from_utf8_lossy(&output.stderr);
                state.add_event(format!("  FAIL  {} — infer_experiment", sample.name));
                warn!(
                    "{}: infer_experiment.py failed (non-fatal): {}",
                    sample.name, stderr
                );
            }
        }

        if crate::is_cancelled() {
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

            match run_cancellable(cmd, "geneBody_coverage.py") {
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

        if crate::is_cancelled() {
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
            let output = cmd.output().map_err(|e| {
                format!("{}: Failed to run {}: {}", sample.name, script.display(), e)
            })?;
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
    for subdir in ["_STARgenome", "_STARpass1", "_STARtmp"] {
        let dir = star_dir.join(format!("{sample_name}{subdir}"));
        if dir.exists() {
            let _ = fs::remove_dir_all(&dir);
        }
    }
}

// ─── Work queue ──────────────────────────────────────────────────────────────

pub(crate) fn run_work_queue<T, F>(
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
                if crate::is_cancelled() {
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

// ─── Genome index detection & generation ─────────────────────────────────────

/// Check if a STAR genome index exists in the given directory.
/// Looks for the `Genome`, `SA`, and `SAindex` files that STAR generates.
pub(crate) fn genome_index_exists(genome_dir: &Path) -> bool {
    genome_dir.join("Genome").exists()
        && genome_dir.join("SA").exists()
        && genome_dir.join("SAindex").exists()
}

/// Find a FASTA reference file (.fa, .fna, .fasta) in the given directory.
/// Returns the first match found.
pub(crate) fn find_fasta_in_dir(dir: &Path) -> Option<PathBuf> {
    let Ok(entries) = fs::read_dir(dir) else {
        return None;
    };
    for entry in entries.filter_map(Result::ok) {
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
            match ext.to_lowercase().as_str() {
                "fa" | "fna" | "fasta" => return Some(path),
                "gz" => {
                    // Check for .fa.gz, .fna.gz, .fasta.gz
                    let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("");
                    if stem.ends_with(".fa") || stem.ends_with(".fna") || stem.ends_with(".fasta") {
                        return Some(path);
                    }
                }
                _ => {}
            }
        }
    }
    None
}

/// Generate a STAR genome index from a FASTA file + GTF annotation.
///
/// Runs: STAR --runMode genomeGenerate --genomeDir <dir> --genomeFastaFiles <fasta>
///       --sjdbGTFfile <gtf> --runThreadN <threads>
///
/// This is a long-running operation (30-90 min for human genome).
pub(crate) fn generate_genome_index(
    config: &Config,
    fasta_path: &Path,
) -> Result<(), String> {
    let star_bin = config.star_env.join("bin/STAR");
    let genome_dir = &config.genome_dir;

    // Ensure genome_dir exists
    fs::create_dir_all(genome_dir).map_err(|e| {
        format!(
            "Cannot create genome index directory {}: {}",
            genome_dir.display(),
            e
        )
    })?;

    // Use all available threads for indexing (it's a single-shot operation)
    let index_threads = config.parallel_jobs * config.threads_per_sample;

    let mut cmd = Command::new(&star_bin);
    cmd.args([
        "--runMode",
        "genomeGenerate",
        "--genomeDir",
        genome_dir.to_str().unwrap(),
        "--genomeFastaFiles",
        fasta_path.to_str().unwrap(),
        "--runThreadN",
        &index_threads.to_string(),
    ]);

    // Add GTF for splice-junction-aware indexing if available
    if config.gtf.exists() {
        cmd.args(["--sjdbGTFfile", config.gtf.to_str().unwrap()]);
    }

    // Log to file in genome_dir
    let log_path = genome_dir.join("genomeGenerate.log");
    let log_file = File::create(&log_path).map_err(|e| {
        format!("Cannot create log file {}: {}", log_path.display(), e)
    })?;
    let log_file2 = log_file.try_clone().unwrap();
    cmd.stdout(Stdio::from(log_file))
        .stderr(Stdio::from(log_file2));

    match run_cancellable(cmd, "STAR genomeGenerate") {
        Ok(true) => {
            // Verify index was created
            if genome_index_exists(genome_dir) {
                Ok(())
            } else {
                Err(format!(
                    "STAR genomeGenerate completed but index files not found in {}\n  \
                     Check log: {}",
                    genome_dir.display(),
                    log_path.display()
                ))
            }
        }
        Ok(false) => Err(format!(
            "STAR genomeGenerate failed.\n  \
             Check log: {}\n  \
             Common causes: insufficient RAM, disk space, or invalid FASTA/GTF files.",
            log_path.display()
        )),
        Err(e) => Err(e),
    }
}

// ─── Environment validation ──────────────────────────────────────────────────

pub(crate) fn validate_environment(config: &Config) -> Result<(), String> {
    let star_bin = config.star_env.join("bin/STAR");
    if !star_bin.exists() {
        return Err(format!(
            "STAR binary not found: {}\n  \
             Hint: Is --star-env correct? Try: --star-env auto\n  \
             Or install: mamba create -n star -c bioconda star=2.7.11b samtools",
            star_bin.display()
        ));
    }

    if !config.samtools.exists() {
        return Err(format!(
            "samtools not found: {}\n  \
             Hint: Specify with --samtools /path/to/samtools\n  \
             Or ensure samtools is installed in the STAR conda env",
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
                "RSeQC script not found: {}\n  \
                 Hint: Is --rseqc-env correct? Try: --rseqc-env auto\n  \
                 Or install: mamba create -n rseqc -c bioconda rseqc python",
                path.display()
            ));
        }
    }

    if !config.genome_dir.exists() {
        return Err(format!(
            "STAR genome dir not found: {}\n  \
             Hint: Specify with --genome-dir /path/to/star/index\n  \
             Or place a reference FASTA (.fa/.fna/.fasta) in the directory for auto-indexing.",
            config.genome_dir.display()
        ));
    }

    // Genome index check is deferred to ensure_genome_index() in main —
    // if a FASTA is present, we can auto-generate the index.

    if !config.gtf.exists() {
        return Err(format!(
            "GTF not found: {}\n  \
             Hint: Specify with --gtf /path/to/annotation.gtf",
            config.gtf.display()
        ));
    }

    if !config.fastq_dir.exists() {
        return Err(format!(
            "FASTQ directory not found: {}",
            config.fastq_dir.display()
        ));
    }

    // Check samtools version (warn only)
    if let Ok(output) = Command::new(&config.samtools).arg("--version").output() {
        if output.status.success() {
            let stdout = String::from_utf8_lossy(&output.stdout);
            if let Some(first_line) = stdout.lines().next() {
                // First line is like "samtools 1.21"
                let version_str = first_line.split_whitespace().nth(1).unwrap_or("");
                let parts: Vec<&str> = version_str.split('.').collect();
                let major: u32 = parts.first().and_then(|s| s.parse().ok()).unwrap_or(0);
                let minor: u32 = parts.get(1).and_then(|s| s.parse().ok()).unwrap_or(0);
                if major < 1 || (major == 1 && minor < 15) {
                    warn!(
                        "samtools version {} detected; version 1.15+ is recommended. \
                         Some features may not work correctly.",
                        version_str
                    );
                }
            }
        }
    }

    Ok(())
}

/// Ensure STAR genome index exists. If not, detect a FASTA file and generate it.
/// Returns Ok(()) if the index is ready, Err if it cannot be created.
pub(crate) fn ensure_genome_index(config: &Config) -> Result<(), String> {
    if genome_index_exists(&config.genome_dir) {
        return Ok(());
    }

    // No index found — look for a FASTA file to generate from
    let fasta_path = find_fasta_in_dir(&config.genome_dir).ok_or_else(|| {
        format!(
            "STAR genome index not found in {} and no reference FASTA \
             (.fa, .fna, .fasta) detected for auto-indexing.\n  \
             Hint: Place your reference genome FASTA in --genome-dir, or run:\n  \
             STAR --runMode genomeGenerate --genomeDir {} --genomeFastaFiles <ref.fa> \
             --sjdbGTFfile <annotation.gtf>",
            config.genome_dir.display(),
            config.genome_dir.display()
        )
    })?;

    info!(
        "STAR genome index not found. Detected FASTA: {}",
        fasta_path.display()
    );
    info!(
        "Generating genome index in {} (this may take 30-90 minutes for a human genome)...",
        config.genome_dir.display()
    );

    generate_genome_index(config, &fasta_path)?;

    info!("Genome index generated successfully.");
    Ok(())
}

// ─── Summary files ───────────────────────────────────────────────────────────

pub(crate) fn write_summary_files(output_dir: &Path, samples: &[Sample]) {
    let star_dir = output_dir.join("star");
    let qc_dir = output_dir.join("qc");

    #[derive(serde::Serialize)]
    struct SummaryRow {
        sample: String,
        sha256: String,
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

    // TSV
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
