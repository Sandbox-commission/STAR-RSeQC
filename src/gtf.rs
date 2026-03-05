use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use log::warn;

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

pub(crate) fn gtf_to_bed12(gtf_path: &Path, bed_path: &Path) -> Result<usize, String> {
    let gtf_file = File::open(gtf_path)
        .map_err(|e| format!("Cannot open GTF {}: {}", gtf_path.display(), e))?;
    let reader = BufReader::new(gtf_file);

    // Collect exons per transcript: (chrom, strand, Vec<(start, end)>)
    let mut transcripts: HashMap<String, (String, String, Vec<(u64, u64)>)> = HashMap::new();
    let mut skipped_exons: usize = 0;

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
            _ => {
                skipped_exons += 1;
                continue;
            }
        };

        transcripts
            .entry(transcript_id)
            .or_insert_with(|| (chrom.to_string(), strand.to_string(), Vec::new()))
            .2
            .push((start, end));
    }

    if skipped_exons > 0 {
        warn!(
            "Skipped {} exon lines without transcript_id in {}",
            skipped_exons,
            gtf_path.display()
        );
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
        return Err("GTF->BED12 produced zero transcripts".to_string());
    }
    Ok(count)
}
