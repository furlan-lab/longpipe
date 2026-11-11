use anyhow::{anyhow, bail, Context, Result};
use clap::Parser;
use rayon::prelude::*;
use serde_json::{Map, Value};
use std::collections::HashMap;
use std::fs;
use std::path::{Component, Path, PathBuf};

/// Materialize files listed in a PacBio MiniWDL output manifest.
#[derive(Parser, Debug)]
#[command(name = "copy_pbWDL", version, about)]
struct Args {
    /// Input JSON result or flat output manifest.
    #[arg(short, long)]
    input: PathBuf,

    /// Directory in which to materialize output files.
    #[arg(short, long)]
    output: PathBuf,

    /// Output manifest (defaults to OUTPUT/outputs.json).
    #[arg(short = 'j', long)]
    json: Option<PathBuf>,

    /// Print planned operations without changing files.
    #[arg(short, long)]
    dry_run: bool,

    /// Maximum number of concurrent file operations.
    #[arg(short = 't', long, value_parser = parse_thread_count)]
    threads: Option<usize>,

    /// Move the haplotagged BAM instead of copying it.
    #[arg(short = 'm', long)]
    move_haplotagged_bam: bool,
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

#[derive(Clone, Debug)]
struct FileOperation {
    key: String,
    source: PathBuf,
    destination: PathBuf,
    move_source: bool,
}

pub fn run() -> Result<()> {
    execute(Args::parse())
}

fn execute(args: Args) -> Result<()> {
    let output_directory = absolute_lexical_path(&args.output)?;
    let json_output = absolute_lexical_path(
        args.json
            .as_deref()
            .unwrap_or_else(|| Path::new("outputs.json")),
    )?;
    let json_output = if args.json.is_some() {
        json_output
    } else {
        output_directory.join("outputs.json")
    };

    let json_content = fs::read_to_string(&args.input)
        .with_context(|| format!("could not read {}", args.input.display()))?;
    let json_data: Value = serde_json::from_str(&json_content)
        .with_context(|| format!("could not parse {}", args.input.display()))?;
    let outputs = workflow_outputs(json_data)?;

    let mut operations = Vec::new();
    let mut destinations = HashMap::new();
    let mut rewritten = Map::new();
    for (key, value) in outputs {
        let value = rewrite_value(
            &key,
            value,
            &output_directory,
            args.move_haplotagged_bam,
            &mut operations,
            &mut destinations,
        )?;
        rewritten.insert(key, value);
    }
    let rewritten = Value::Object(rewritten);

    if args.dry_run {
        println!("Dry run - planned file operations:");
        for operation in &operations {
            let verb = if operation.move_source {
                "Move"
            } else {
                "Copy"
            };
            println!(
                "  {verb} {} to {} ({})",
                operation.source.display(),
                operation.destination.display(),
                operation.key,
            );
        }
        println!("\nWould write manifest to {}", json_output.display());
        println!("{}", serde_json::to_string_pretty(&rewritten)?);
        return Ok(());
    }

    fs::create_dir_all(&output_directory).with_context(|| {
        format!(
            "could not create output directory {}",
            output_directory.display()
        )
    })?;

    let thread_count = args
        .threads
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, usize::from));
    println!(
        "Materializing {} files with {} threads",
        operations.len(),
        thread_count
    );
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build()
        .context("could not create copy thread pool")?;
    pool.install(|| {
        operations
            .par_iter()
            .enumerate()
            .map(|(index, operation)| perform_operation(operation, index))
            .collect::<Result<Vec<_>>>()
    })?;

    write_json_atomically(&json_output, &rewritten)?;
    println!("Wrote manifest to {}", json_output.display());
    Ok(())
}

fn workflow_outputs(json_data: Value) -> Result<Map<String, Value>> {
    let mut root = match json_data {
        Value::Object(map) => map,
        _ => bail!("JSON root must be an object"),
    };

    let is_miniwdl_result =
        root.get("dir").is_some_and(Value::is_string) && root.contains_key("outputs");
    if !is_miniwdl_result {
        return Ok(root);
    }

    match root.remove("outputs") {
        Some(Value::Object(outputs)) => Ok(outputs),
        Some(_) => bail!("MiniWDL result field 'outputs' must be an object"),
        None => unreachable!("outputs was checked above"),
    }
}

fn rewrite_value(
    key: &str,
    value: Value,
    output_directory: &Path,
    move_haplotagged_bam: bool,
    operations: &mut Vec<FileOperation>,
    destinations: &mut HashMap<PathBuf, PathBuf>,
) -> Result<Value> {
    match value {
        Value::String(string_value) if is_likely_file_path(&string_value) => {
            let source = PathBuf::from(&string_value);
            if !source.exists() {
                eprintln!(
                    "Warning: path does not exist for {key}: {}",
                    source.display()
                );
                return Ok(Value::String(string_value));
            }
            let source = fs::canonicalize(&source).with_context(|| {
                format!("could not resolve path for {key}: {}", source.display())
            })?;
            let file_name = source
                .file_name()
                .ok_or_else(|| anyhow!("source has no filename: {}", source.display()))?;
            let destination = output_directory.join(file_name);

            match destinations.get(&destination) {
                Some(existing_source) if existing_source != &source => {
                    bail!(
                        "output filename collision: {} and {} both map to {}",
                        existing_source.display(),
                        source.display(),
                        destination.display()
                    );
                }
                Some(_) => {}
                None => {
                    destinations.insert(destination.clone(), source.clone());
                    operations.push(FileOperation {
                        key: key.to_string(),
                        move_source: move_haplotagged_bam && is_haplotagged_bam(key, &source),
                        source,
                        destination: destination.clone(),
                    });
                }
            }
            Ok(Value::String(destination.to_string_lossy().into_owned()))
        }
        Value::Array(values) => values
            .into_iter()
            .map(|value| {
                rewrite_value(
                    key,
                    value,
                    output_directory,
                    move_haplotagged_bam,
                    operations,
                    destinations,
                )
            })
            .collect::<Result<Vec<_>>>()
            .map(Value::Array),
        other => Ok(other),
    }
}

fn is_likely_file_path(value: &str) -> bool {
    value.starts_with('/')
        || value.contains('/')
        || Path::new(value)
            .file_name()
            .and_then(|name| name.to_str())
            .is_some_and(|name| name.contains('.'))
}

fn is_haplotagged_bam(key: &str, source: &Path) -> bool {
    let filename = source
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or_default();
    (filename.contains("haplotagged") && filename.ends_with(".bam"))
        || key.contains("haplotagged_bam")
}

fn perform_operation(operation: &FileOperation, index: usize) -> Result<()> {
    let parent = operation.destination.parent().ok_or_else(|| {
        anyhow!(
            "destination has no parent: {}",
            operation.destination.display()
        )
    })?;
    fs::create_dir_all(parent).with_context(|| format!("could not create {}", parent.display()))?;

    if operation.move_source {
        match fs::rename(&operation.source, &operation.destination) {
            Ok(()) => return Ok(()),
            Err(error) if is_cross_device_error(&error) => {
                copy_atomically(operation, index)?;
                fs::remove_file(&operation.source).with_context(|| {
                    format!(
                        "copied but could not remove moved source {}",
                        operation.source.display()
                    )
                })?;
                return Ok(());
            }
            Err(error) => {
                return Err(error).with_context(|| {
                    format!(
                        "could not move {} to {} for {}",
                        operation.source.display(),
                        operation.destination.display(),
                        operation.key
                    )
                });
            }
        }
    }

    copy_atomically(operation, index)
}

fn copy_atomically(operation: &FileOperation, index: usize) -> Result<()> {
    let filename = operation
        .destination
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("output");
    let temporary = operation.destination.with_file_name(format!(
        ".{filename}.longpipe-copy-{}-{index}",
        std::process::id()
    ));
    let result = (|| {
        fs::copy(&operation.source, &temporary).with_context(|| {
            format!(
                "could not copy {} for {}",
                operation.source.display(),
                operation.key
            )
        })?;
        fs::rename(&temporary, &operation.destination).with_context(|| {
            format!(
                "could not publish copied file {}",
                operation.destination.display()
            )
        })?;
        Ok(())
    })();
    if result.is_err() {
        let _ = fs::remove_file(&temporary);
    }
    result
}

#[cfg(unix)]
fn is_cross_device_error(error: &std::io::Error) -> bool {
    error.raw_os_error() == Some(18)
}

#[cfg(not(unix))]
fn is_cross_device_error(_error: &std::io::Error) -> bool {
    false
}

fn write_json_atomically(path: &Path, value: &Value) -> Result<()> {
    let parent = path
        .parent()
        .ok_or_else(|| anyhow!("manifest has no parent: {}", path.display()))?;
    fs::create_dir_all(parent).with_context(|| format!("could not create {}", parent.display()))?;
    let filename = path
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("outputs.json");
    let temporary = path.with_file_name(format!(
        ".{filename}.longpipe-manifest-{}",
        std::process::id()
    ));
    let result = (|| {
        let mut bytes = serde_json::to_vec_pretty(value)?;
        bytes.push(b'\n');
        fs::write(&temporary, bytes)
            .with_context(|| format!("could not write {}", temporary.display()))?;
        fs::rename(&temporary, path)
            .with_context(|| format!("could not publish {}", path.display()))?;
        Ok(())
    })();
    if result.is_err() {
        let _ = fs::remove_file(&temporary);
    }
    result
}

fn absolute_lexical_path(path: &Path) -> Result<PathBuf> {
    let absolute = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir()
            .context("could not determine current directory")?
            .join(path)
    };
    let mut normalized = PathBuf::new();
    for component in absolute.components() {
        match component {
            Component::CurDir => {}
            Component::ParentDir => {
                normalized.pop();
            }
            other => normalized.push(other.as_os_str()),
        }
    }
    Ok(normalized)
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn temporary_directory(name: &str) -> PathBuf {
        let directory = std::env::temp_dir().join(format!(
            "longpipe-{name}-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ));
        if directory.exists() {
            fs::remove_dir_all(&directory).unwrap();
        }
        fs::create_dir_all(&directory).unwrap();
        directory
    }

    #[test]
    fn extracts_outputs_from_miniwdl_result() {
        let outputs = workflow_outputs(json!({
            "dir": "/tmp/run/workflow",
            "outputs": {
                "workflow.bam": "/tmp/run/workflow/sample.bam",
                "workflow.version": "1.0"
            }
        }))
        .unwrap();
        assert_eq!(outputs.len(), 2);
        assert!(!outputs.contains_key("dir"));
    }

    #[test]
    fn preserves_flat_manifest() {
        let outputs = workflow_outputs(json!({
            "workflow.bam": "/tmp/run/workflow/sample.bam",
            "workflow.version": "1.0"
        }))
        .unwrap();
        assert_eq!(outputs.len(), 2);
        assert_eq!(outputs.get("workflow.version"), Some(&json!("1.0")));
    }

    #[test]
    fn rejects_malformed_miniwdl_result() {
        let error = workflow_outputs(json!({
            "dir": "/tmp/run/workflow",
            "outputs": []
        }))
        .unwrap_err();
        assert_eq!(
            error.to_string(),
            "MiniWDL result field 'outputs' must be an object"
        );
    }

    #[test]
    fn detects_paths_without_misclassifying_versions() {
        assert!(is_likely_file_path("sample.bam"));
        assert!(is_likely_file_path("results/sample"));
        assert!(!is_likely_file_path("v2"));
    }

    #[test]
    fn materializes_miniwdl_files_and_moves_the_haplotagged_bam() {
        let root = temporary_directory("copy-pbwdl");
        let source_directory = root.join("source");
        let output_directory = root.join("output");
        fs::create_dir(&source_directory).unwrap();
        let bam = source_directory.join("sample.GRCh38.haplotagged.bam");
        let report = source_directory.join("sample.report.txt");
        fs::write(&bam, b"bam").unwrap();
        fs::write(&report, b"report").unwrap();
        let manifest = root.join("miniwdl.json");
        fs::write(
            &manifest,
            serde_json::to_vec_pretty(&json!({
                "dir": source_directory,
                "outputs": {
                    "humanwgs_singleton.merged_haplotagged_bam": bam,
                    "humanwgs_singleton.reports": [[report]],
                    "humanwgs_singleton.workflow_version": "v2.0.7"
                }
            }))
            .unwrap(),
        )
        .unwrap();

        execute(Args {
            input: manifest,
            output: output_directory.clone(),
            json: None,
            dry_run: false,
            threads: Some(2),
            move_haplotagged_bam: true,
        })
        .unwrap();

        assert!(!bam.exists());
        assert_eq!(
            fs::read(output_directory.join("sample.GRCh38.haplotagged.bam")).unwrap(),
            b"bam"
        );
        assert!(report.exists());
        assert_eq!(
            fs::read(output_directory.join("sample.report.txt")).unwrap(),
            b"report"
        );
        let rewritten: Value =
            serde_json::from_slice(&fs::read(output_directory.join("outputs.json")).unwrap())
                .unwrap();
        assert_eq!(rewritten["humanwgs_singleton.workflow_version"], "v2.0.7");
        assert_eq!(
            rewritten["humanwgs_singleton.merged_haplotagged_bam"],
            output_directory
                .join("sample.GRCh38.haplotagged.bam")
                .to_string_lossy()
                .as_ref()
        );
        assert!(rewritten.get("dir").is_none());

        fs::remove_dir_all(root).unwrap();
    }
}
