use anyhow::{anyhow, bail, Context, Result};
use clap::Parser;
use flate2::read::MultiGzDecoder;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

/// Convert a FIRE pileup to a compact four-column bedGraph.
#[derive(Parser, Debug)]
#[command(name = "shrink_fire_bed", version, about)]
struct Args {
    /// Tabix-indexed, bgzipped FIRE pileup BED.
    input: PathBuf,

    /// Destination bedGraph.
    output: PathBuf,

    /// Tab-delimited chromosome sizes file.
    #[arg(short, long)]
    chrom_sizes: PathBuf,

    /// Accepted for compatibility; conversion uses one streaming pass.
    #[arg(short, long, value_parser = parse_thread_count)]
    threads: Option<usize>,
}

fn parse_thread_count(value: &str) -> std::result::Result<usize, String> {
    let threads = value
        .parse::<usize>()
        .map_err(|_| "thread count must be a positive integer".to_string())?;
    if threads == 0 {
        return Err("thread count must be at least 1".to_string());
    }
    Ok(threads)
}

#[derive(Debug, Clone, PartialEq)]
struct BedGraphEntry {
    chrom: String,
    start: u64,
    end: u64,
    value: f64,
}

pub fn run() -> Result<()> {
    let args = Args::parse();
    let chromosome_sizes = read_chromosome_sizes(&args.chrom_sizes)?;
    let reader = open_input(&args.input)?;
    let parent = args
        .output
        .parent()
        .ok_or_else(|| anyhow!("output has no parent directory: {}", args.output.display()))?;
    if !parent.as_os_str().is_empty() {
        fs::create_dir_all(parent)
            .with_context(|| format!("could not create {}", parent.display()))?;
    }

    let filename = args
        .output
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("fire.bedgraph");
    let temporary = args.output.with_file_name(format!(
        ".{filename}.longpipe-shrink-{}",
        std::process::id()
    ));
    let result: Result<usize> = (|| {
        let output = File::create(&temporary)
            .with_context(|| format!("could not create {}", temporary.display()))?;
        let mut writer = BufWriter::new(output);
        let records = shrink(reader, &mut writer, &chromosome_sizes)?;
        writer.flush().context("could not flush output bedGraph")?;
        drop(writer);
        fs::rename(&temporary, &args.output)
            .with_context(|| format!("could not publish {}", args.output.display()))?;
        Ok(records)
    })();
    if result.is_err() {
        let _ = fs::remove_file(&temporary);
    }
    let records = result?;
    let requested_threads = args.threads.unwrap_or(1);
    eprintln!(
        "Wrote {records} bedGraph records to {} (streaming; requested threads: {requested_threads})",
        args.output.display()
    );
    Ok(())
}

fn open_input(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path)
        .with_context(|| format!("could not open FIRE pileup {}", path.display()))?;
    let compressed = matches!(
        path.extension().and_then(|extension| extension.to_str()),
        Some("gz" | "bgz" | "bgzf")
    );
    if compressed {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn read_chromosome_sizes(path: &Path) -> Result<HashMap<String, u64>> {
    let file = File::open(path)
        .with_context(|| format!("could not open chromosome sizes {}", path.display()))?;
    let mut sizes = HashMap::new();
    for (index, line) in BufReader::new(file).lines().enumerate() {
        let line_number = index + 1;
        let line =
            line.with_context(|| format!("could not read {}:{line_number}", path.display()))?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let mut fields = line.split_whitespace();
        let chromosome = fields
            .next()
            .ok_or_else(|| anyhow!("{}:{line_number}: missing chromosome", path.display()))?;
        let length = fields
            .next()
            .ok_or_else(|| {
                anyhow!(
                    "{}:{line_number}: missing length for {chromosome}",
                    path.display()
                )
            })?
            .parse::<u64>()
            .with_context(|| {
                format!(
                    "{}:{line_number}: invalid length for {chromosome}",
                    path.display()
                )
            })?;
        if length == 0 {
            bail!(
                "{}:{line_number}: chromosome {chromosome} has zero length",
                path.display()
            );
        }
        if sizes.insert(chromosome.to_string(), length).is_some() {
            bail!(
                "{}:{line_number}: duplicate chromosome {chromosome}",
                path.display()
            );
        }
    }
    if sizes.is_empty() {
        bail!("chromosome sizes file is empty: {}", path.display());
    }
    Ok(sizes)
}

fn shrink<R: BufRead, W: Write>(
    reader: R,
    writer: &mut W,
    chromosome_sizes: &HashMap<String, u64>,
) -> Result<usize> {
    let mut current: Option<BedGraphEntry> = None;
    let mut completed_chromosomes = HashSet::new();
    let mut active_chromosome: Option<String> = None;
    let mut previous_start = 0;
    let mut written = 0;

    for (index, line) in reader.lines().enumerate() {
        let line_number = index + 1;
        let line = line.with_context(|| format!("could not read pileup line {line_number}"))?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 6 {
            bail!(
                "pileup line {line_number} has {} columns; expected at least 6",
                fields.len()
            );
        }
        let chromosome = fields[0];
        let chromosome_length = chromosome_sizes.get(chromosome).ok_or_else(|| {
            anyhow!(
                "pileup line {line_number}: chromosome {chromosome} is absent from chromosome sizes"
            )
        })?;
        let start = fields[1]
            .parse::<u64>()
            .with_context(|| format!("pileup line {line_number}: invalid start"))?;
        let end = fields[2]
            .parse::<u64>()
            .with_context(|| format!("pileup line {line_number}: invalid end"))?;
        let value = fields[5]
            .parse::<f64>()
            .with_context(|| format!("pileup line {line_number}: invalid FIRE score"))?;
        if start >= end {
            bail!("pileup line {line_number}: interval start {start} must be less than end {end}");
        }
        if end > *chromosome_length {
            bail!(
                "pileup line {line_number}: interval ends at {end}, beyond {chromosome} length {chromosome_length}"
            );
        }
        if !value.is_finite() {
            bail!("pileup line {line_number}: FIRE score must be finite");
        }

        match active_chromosome.as_deref() {
            Some(active) if active == chromosome => {
                if start < previous_start {
                    bail!("pileup line {line_number}: {chromosome} positions are not sorted");
                }
            }
            Some(active) => {
                completed_chromosomes.insert(active.to_string());
                if completed_chromosomes.contains(chromosome) {
                    bail!("pileup line {line_number}: chromosome {chromosome} is not contiguous");
                }
                active_chromosome = Some(chromosome.to_string());
            }
            None => active_chromosome = Some(chromosome.to_string()),
        }
        previous_start = start;

        let next = BedGraphEntry {
            chrom: chromosome.to_string(),
            start,
            end,
            value,
        };
        if let Some(mut entry) = current.take() {
            if entry.chrom == next.chrom
                && entry.end == next.start
                && entry.value == -1.0
                && next.value == -1.0
            {
                entry.end = next.end;
                current = Some(entry);
            } else {
                write_entry(writer, &entry)?;
                written += 1;
                current = Some(next);
            }
        } else {
            current = Some(next);
        }
    }

    if let Some(entry) = current {
        write_entry(writer, &entry)?;
        written += 1;
    }
    Ok(written)
}

fn write_entry<W: Write>(writer: &mut W, entry: &BedGraphEntry) -> Result<()> {
    let value = if entry.value == -1.0 {
        0.0
    } else {
        entry.value
    };
    writeln!(
        writer,
        "{}\t{}\t{}\t{}",
        entry.chrom, entry.start, entry.end, value
    )
    .context("could not write output bedGraph")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn sizes() -> HashMap<String, u64> {
        HashMap::from([("chr1".to_string(), 100), ("chr2".to_string(), 50)])
    }

    #[test]
    fn merges_adjacent_sentinel_intervals_and_converts_them_to_zero() {
        let input = concat!(
            "#chrom\tstart\tend\tFDR\tis_local_max\tscore\n",
            "chr1\t0\t10\t1\t0\t-1.0\n",
            "chr1\t10\t20\t1\t0\t-1\n",
            "chr1\t20\t30\t0.1\t1\t2.5\n",
            "chr2\t0\t5\t1\t0\t-1\n",
        );
        let mut output = Vec::new();
        let records = shrink(Cursor::new(input), &mut output, &sizes()).unwrap();
        assert_eq!(records, 3);
        assert_eq!(
            String::from_utf8(output).unwrap(),
            concat!("chr1\t0\t20\t0\n", "chr1\t20\t30\t2.5\n", "chr2\t0\t5\t0\n",)
        );
    }

    #[test]
    fn rejects_unknown_chromosomes() {
        let input = "chr3\t0\t10\t1\t0\t-1\n";
        let error = shrink(Cursor::new(input), &mut Vec::new(), &sizes()).unwrap_err();
        assert!(error.to_string().contains("absent from chromosome sizes"));
    }

    #[test]
    fn rejects_noncontiguous_chromosome_blocks() {
        let input = concat!(
            "chr1\t0\t10\t1\t0\t-1\n",
            "chr2\t0\t10\t1\t0\t-1\n",
            "chr1\t10\t20\t1\t0\t-1\n",
        );
        let error = shrink(Cursor::new(input), &mut Vec::new(), &sizes()).unwrap_err();
        assert!(error.to_string().contains("is not contiguous"));
    }
}
