use std::fs::{self, File};
use std::path::{Path, PathBuf};

pub(crate) fn checkpoint_dir(output_dir: &Path) -> PathBuf {
    output_dir.join(".checkpoints")
}

pub(crate) fn sha256_file(path: &Path) -> Result<Vec<u8>, String> {
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

pub(crate) fn sha256_file_list(files: &[PathBuf]) -> String {
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

pub(crate) struct SampleDigests {
    pub star: String,
    pub rseqc: String,
}

pub(crate) fn sha256_outputs(output_dir: &Path, sample_name: &str) -> SampleDigests {
    let star_dir = output_dir.join("star");
    let qc_dir = output_dir.join("qc");

    let star_files: Vec<PathBuf> = vec![
        star_dir.join(format!("{sample_name}_Log.out")),
        star_dir.join(format!("{sample_name}_Log.progress.out")),
        star_dir.join(format!("{sample_name}_Log.final.out")),
        star_dir.join(format!("{sample_name}_Aligned.sortedByCoord.out.bam")),
        star_dir.join(format!(
            "{sample_name}_Aligned.sortedByCoord.out.bam.bai"
        )),
        star_dir.join(format!(
            "{sample_name}_Aligned.toTranscriptome.out.bam"
        )),
        star_dir.join(format!("{sample_name}_ReadsPerGene.out.tab")),
        star_dir.join(format!("{sample_name}_SJ.out.tab")),
        star_dir.join(format!("{sample_name}_Chimeric.out.junction")),
        star_dir.join(format!("{sample_name}_Chimeric.out.sam")),
    ];

    let rseqc_files: Vec<PathBuf> = vec![
        qc_dir.join(format!("{sample_name}.strand.txt")),
        qc_dir.join(format!("{sample_name}.geneBodyCoverage.txt")),
        qc_dir.join(format!("{sample_name}.geneBodyCoverage.r")),
        qc_dir.join(format!(
            "{sample_name}.geneBodyCoverage.curves.pdf"
        )),
        qc_dir.join(format!(
            "{sample_name}.geneBodyCoverage.heatMap.pdf"
        )),
        qc_dir.join(format!("{sample_name}.read_distribution.txt")),
    ];

    SampleDigests {
        star: sha256_file_list(&star_files),
        rseqc: sha256_file_list(&rseqc_files),
    }
}

pub(crate) fn write_checkpoint(output_dir: &Path, name: &str, digests: &SampleDigests) {
    let dir = checkpoint_dir(output_dir);
    let _ = fs::create_dir_all(&dir);
    let path = dir.join(format!("{name}.sha256"));
    let _ = fs::write(
        &path,
        format!("star:{}\nrseqc:{}\n", digests.star, digests.rseqc),
    );
}

pub(crate) fn remove_checkpoint(output_dir: &Path, name: &str) {
    let path = checkpoint_dir(output_dir).join(format!("{name}.sha256"));
    let _ = fs::remove_file(&path);
}

pub(crate) fn parse_checkpoint(content: &str) -> Option<(String, String)> {
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

pub(crate) enum ResumeStatus {
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

pub(crate) fn check_resume(output_dir: &Path, sample_name: &str) -> ResumeStatus {
    let ckpt = checkpoint_dir(output_dir).join(format!("{sample_name}.sha256"));
    let content = match fs::read_to_string(&ckpt) {
        Ok(s) => s,
        Err(_) => return ResumeStatus::NotDone,
    };

    let (old_star, old_rseqc) = match parse_checkpoint(&content) {
        Some(pair) => pair,
        None => return ResumeStatus::NotDone,
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
