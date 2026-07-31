# Longpipe Development Notes

## Repository Purpose

Longpipe is a Snakemake workflow for PacBio HiFi and Fiber-seq samples. The
workflow currently runs:

- optional Fiber-seq nucleosome annotation;
- PacBio HiFi human WGS WDL through MiniWDL;
- FIRE, bigWig generation, and cCRE aggregation for Hia5 samples;
- Sawfish structural-variant calling and SVTopo extraction;
- HiFiHLA from phased contigs or aligned reads.

The Python CLI is implemented in `src/longpipe/cli.py`. Workflow resources
remain at the repository root during development and are installed under
`share/longpipe` by `pyproject.toml`. `Snakefile` requires the execution
context supplied by the CLI.

## CLI Contract

The user creates a project working directory, normally below a configurable
fast-storage root such as:

```text
/fh/fast/furlan_s/project/
```

From a chosen working directory, the user initializes a named run:

```bash
longpipe setup RUNNAME
```

This creates a durable run directory:

```text
<working-directory>/
└── RUNNAME/
    ├── config.yaml
    └── runsheet.csv
```

`setup` copies editable templates into that directory. It must refuse to
overwrite an existing run directory.

After editing the files, the user starts the workflow from the same working
directory:

```bash
longpipe run RUNNAME
```

The command remains attached until Snakemake completes or fails. Input data,
including source BAM files on fast storage, is read in place and is not copied
to scratch.

## Storage Model

The complete Snakemake run takes place under a configurable scratch root, whose
site default is:

```text
/hpc/temp/furlan_s/
```

For a run named `RUNNAME`, its scratch working directory is:

```text
/hpc/temp/furlan_s/RUNNAME/
```

The fast-storage project root and scratch root must be configurable rather than
hard-coded to one user. The CLI should validate that setup and run operations
are being performed below the configured project root.

The CLI derives the workflow name and scratch work directory from `RUNNAME`.
They should not be independently editable in `config.yaml`; specifically,
`name` and `outdir` should be supplied by the CLI.

Relative paths in `config.yaml` and `runsheet.csv` are resolved relative to the
durable run directory containing those files. Absolute paths remain valid.

## Configuration Snapshot

At launch, copy `config.yaml` and `runsheet.csv` from the durable run directory
into the scratch run directory. The scratch copies are the immutable inputs
for that invocation. Edits made to the durable copies after launch must not
alter a running workflow.

Published provenance must include the exact configuration and runsheet used by
the workflow, not potentially newer edits from the durable run directory.

## Published Files

On successful completion, copy final scientific outputs and operational
records from scratch back to the durable run directory. Use this general
layout:

```text
RUNNAME/
├── config.yaml
├── runsheet.csv
├── results/
│   ├── fire/
│   ├── sv/
│   └── hla/
├── logs/
├── benchmarks/
└── run metadata/status
```

FIRE's user-facing outputs must be published below `results/fire/`, rather
than retaining the workflow's current top-level `fibertools/FIRE` paths.
Publish the final FIRE pileups, bigWigs, and cCRE summaries.

`ccre_peaks` is a non-empty YAML mapping of filesystem-safe peak-set names to
BED paths. Run cCRE aggregation for every Fiber-seq sample and every configured
peak set, publishing each set below `results/fire/ccre/<peak-set>/`.

Publish final Sawfish/SVTopo and HiFiHLA outputs below `results/sv/` and
`results/hla/`, respectively. Also publish logs, benchmarks, Snakemake
metadata needed for diagnosis, the configuration snapshot, and a clear run
status record.

Large intermediates remain only in scratch. These include MiniWDL run
directories, assemblies, temporary FASTQs, Sawfish working directories,
Conda environments, and other non-final workflow products.

If the run fails, preserve scratch and publish diagnostics: logs, benchmarks,
the configuration snapshot, and a failure/status record. Do not publish
partial scientific results as though they were completed outputs.

## Resume, Collision, and Cleanup Semantics

- Scratch is retained after both successful and failed runs.
- `longpipe run RUNNAME` resumes an existing scratch run by default.
- `longpipe run RUNNAME --restart` is required to discard existing scratch
  state and start over.
- Restart must never happen implicitly.
- Concurrent runs with the same run name must be rejected, including runs
  launched from different project working directories, because their scratch
  paths would collide.
- Use a reliable run lock and record enough metadata to identify its owning
  invocation.
- Provide `longpipe clean RUNNAME` as an explicit way to remove retained
  scratch data. Cleanup must not remove durable published results.

## Implementation Constraints

- Keep Snakemake as the workflow engine and continue using the SLURM profile.
- Pass the selected configuration to Snakemake explicitly; do not depend on
  the repository's template config at execution time.
- Ensure schema validation applies to the per-run configuration and
  runsheet.
- Do not copy source BAMs or reference data into scratch.
- Publish files only after their producing workflow rules have completed.
- Use atomic or staged publication where practical so interrupted copying
  cannot make an incomplete result set appear successful.
- Preserve enough scratch state for Snakemake to resume without rerunning
  completed work unnecessarily.
