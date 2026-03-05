use glob::glob;
use log::warn;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone)]
pub(crate) struct Sample {
    pub name: String,
    pub r1: PathBuf,
    pub r2: PathBuf,
}

pub(crate) fn discover_samples(fastq_dir: &Path, r1_suffix: &str, r2_suffix: &str) -> Vec<Sample> {
    let pattern = fastq_dir
        .join(format!("*{r1_suffix}.fastq.gz"))
        .to_string_lossy()
        .to_string();

    let suffix_with_ext = format!("{r1_suffix}.fastq.gz");

    let mut samples = Vec::new();
    let mut seen = HashMap::new();

    let entries: Vec<_> = match glob(&pattern) {
        Ok(paths) => paths.filter_map(|e| e.ok()).collect(),
        Err(_) => return samples,
    };

    for r1 in entries {
        let r1_name = r1.file_name().unwrap().to_string_lossy().to_string();

        let sample_name = match r1_name.strip_suffix(&suffix_with_ext) {
            Some(n) => n.to_string(),
            None => continue,
        };

        let r2_name = format!("{sample_name}{r2_suffix}.fastq.gz");
        let r2 = r1.parent().unwrap().join(&r2_name);

        if !r2.exists() {
            warn!("Skipping {} — R2 not found ({})", sample_name, r2.display());
            continue;
        }

        if seen.contains_key(&sample_name) {
            warn!("Duplicate sample name {} — skipping", sample_name);
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
