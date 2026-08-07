# longpipe

`longpipe` runs a Snakemake workflow for PacBio HiFi and Fiber-seq samples. It
includes:

1. optional Fiber-seq nucleosome annotation;
2. PacBio HiFi human WGS WDL through MiniWDL;
3. FIRE, haplotype accessibility tracks and tests, and cCRE aggregation for
   Hia5 samples;
4. Sawfish structural-variant calling and SVTopo extraction;
5. HiFiHLA from phased contigs or aligned reads.

The top-level Snakemake process stays attached on the submit host. MiniWDL is
the only nested workflow driver. FIRE pileup generation and the remaining
rules are submitted directly by Snakemake's SLURM executor.

## Installation

Install the command in an isolated pipx environment using Python 3.12:

```bash
ml Python/3.12.3-GCCcore-13.3.0
pipx install --python python3.12 .
longpipe --version
```

This installs Longpipe, Snakemake, pandas, PyYAML, and the Snakemake SLURM
executor without requiring `conda activate longpipe`. Run `pipx ensurepath`
first if pipx's application directory is not already on `PATH`.

### Python Module Bootstrap

On a module-based HPC system, an existing `pipx` command may point to a
different Python module. If pipx reports a missing `libpython` for an older
Python version, bootstrap pipx using the same Python 3.12 module that will run
Longpipe:

```bash
ml Python/3.12.3-GCCcore-13.3.0

python3 -m venv "$HOME/.local/share/pipx-3.12"
"$HOME/.local/share/pipx-3.12/bin/python" -m pip install --upgrade pip pipx
mkdir -p "$HOME/.local/bin"
ln -sf "$HOME/.local/share/pipx-3.12/bin/pipx" "$HOME/.local/bin/pipx"
"$HOME/.local/share/pipx-3.12/bin/pipx" ensurepath
hash -r

pipx install --python "$(command -v python3)" .
longpipe --version
```

Load `Python/3.12.3-GCCcore-13.3.0` before invoking pipx or Longpipe in future
shells. Both virtual environments reference that module's Python installation.

Conda is still required as a workflow deployment tool. The SLURM profile uses
the environment definitions in `envs/` to create isolated Bioconda
environments for tools such as fibertools-rs, hifiasm, HiFiHLA, samtools,
bedtools, pbtk, Sawfish, and SVTopo. A working `conda` executable must
therefore be available on the submit and compute nodes, but no Longpipe Conda
environment needs to be activated.

MiniWDL remains externally managed because it requires site scheduler
configuration. FIRE's pileup-producing DAG is maintained in Longpipe and uses
the pinned per-rule environments in `envs/fire-core.yaml` and
`envs/fire-python.yaml`; users do not install or run a shared FIRE Pixi
environment. The `copy_pbWDL` and `shrink_fire_bed` executables from
pacbiowdlR remain external. Set their executable and resource paths in each
run's `config.yaml`.

The native FIRE rules reproduce the pileup path from upstream FIRE commit
`e7f2aae1f59bbb74c46e1f27936da380cb074704` (`v0.1.2`). Vendored FIRE code is
covered by `workflow/fire/LICENSE`. `fire_exclude_regions` must list the gap,
blacklist, and segmental-duplication BED inputs used by that release.
Existing run configurations must remove `fire_pixi` and `pixi_executable`,
then add `fire_upstream_commit` and `fire_exclude_regions` as shown in
`config/config.yaml` before testing this branch.

The Furlan Lab MiniWDL deployment files are in `deploy/miniwdl/`. Its wrapper
uses a shared engine and site configuration while keeping call and Singularity
caches below the configured Longpipe scratch root. Set `MINIWDL_CACHE_ROOT` to
use a different persistent cache location.

Sawfish's `annotation_and_common_cnv` and sex-specific expected-copy-number
BED files are reference inputs rather than software dependencies. Keep them in
a shared, versioned reference directory matching the Sawfish version declared
in `envs/sawfish.yaml`.

The WDL checkout, reference map, and tertiary map must come from the same
HiFi-human-WGS-WDL resource release.

## Storage Configuration

Longpipe separates durable project data from temporary workflow state. The
site defaults are:

```text
project root: /fh/fast/furlan_s/project
scratch root: /hpc/temp/furlan_s
```

Override them for another user or installation:

```bash
export LONGPIPE_PROJECT_ROOT=/fh/fast/$USER/project
export LONGPIPE_SCRATCH_ROOT=/hpc/temp/$USER
```

The same values can be passed as `--project-root` and `--scratch-root` on the
relevant commands. Longpipe rejects operations launched outside the configured
project root.

## Create A Run

Create or choose a working directory below the project root, then initialize a
named run:

```bash
cd /fh/fast/furlan_s/project/my-project
longpipe setup 260410_BT-brain
```

This creates:

```text
260410_BT-brain/
├── config.yaml
└── runsheet.csv
```

`setup` refuses to overwrite an existing directory. Edit both files before
starting the run. The CLI derives the run name and scratch location, so
`config.yaml` does not contain `name` or `outdir`.

Relative paths in `config.yaml` and `runsheet.csv` resolve from the durable run
directory. Absolute paths are also supported. Source BAMs and reference data
are read in place rather than staged into scratch.

Configure one or more named cCRE peak sets as a YAML mapping:

```yaml
ccre_peaks:
  encode: refs/encode_cCREs.bed.gz
  brain: /fh/fast/shared/brain_cCREs.bed.gz
```

Each peak set is mapped against every Fiber-seq sample. Names may contain
letters, digits, `.`, `_`, and `-`, and are used as output directory names
below `samples/<sample>/fire/ccre/`.

FIRE runs only for Fiber-seq samples whose PacBio WDL `outputs.json` reports
the PacBio WDL mean-depth output at or above `fire_min_mean_depth`. Longpipe
prefers the current `humanwgs_singleton.stat_depth_mean` key and also accepts
the legacy `humanwgs_singleton.stat_mean_depth` key. The threshold defaults
to `8.0`. This is a mean-depth prefilter; FIRE's own coverage checks
still apply. Decisions and skipped samples are recorded in
`run-metadata/fire-eligibility.tsv`.

FIRE's filtering and coverage defaults are explicit in `config.yaml`. Keep
`fire_upstream_commit` and `fire_version` pinned unless the native rules and
parity tests are updated together. To compare an old nested-FIRE result with a
native result while testing this branch, compare decompressed content:

```bash
python workflow/scripts/compare_fire_pileups.py \
  /path/to/legacy-pileup.bed.gz \
  /path/to/native-pileup.bed.gz
```

The final pileup includes haplotype-specific coverage and FIRE coverage from
`ft pileup --haps`. Longpipe converts those columns into H1 and H2 percent
accessible bigWigs plus an H1-minus-H2 percentage-point bigWig. It also tests
confident FIRE local maxima for a haplotype accessibility difference using a
two-sided Fisher exact test on accessible versus inaccessible fibers and
Benjamini-Hochberg correction across tested loci. The thresholds are controlled
by `fire_haplotype_min_coverage`, `fire_haplotype_max_fdr`, and
`fire_haplotype_qvalue`.

This is a phased H1-versus-H2 analysis, not a cCRE aggregation. H1 and H2 are
phase labels and are not automatically maternal and paternal. Their direction
can switch between unoriented phase blocks, and unphased fibers do not
contribute to either haplotype-specific denominator.

The runsheet requires:

| Column | Meaning |
| --- | --- |
| `core_name` | Core identifier and fallback sample label |
| `lab_name` | Optional preferred sample label |
| `run` | PacBio run directory below `rawdir` |
| `cell` | SMRT Cell directory |
| `barcode` | Barcode used in the demultiplexed BAM name |
| `hifi_bam` | Optional explicit BAM path; overrides run/cell/barcode discovery |
| `sex` | `M`, `MALE`, `F`, `FEMALE`, or blank for `default_sex` |
| `hia5` | Explicit yes/no value controlling the Fiber-seq branch |

Labels must be unique and contain only letters, digits, `.`, `_`, and `-`.
Configuration and runsheet schemas are checked before jobs launch.

### HLA Modes

- `hla_mode: contigs` preserves broad HiFiHLA locus coverage but performs an
  expensive whole-genome hifiasm assembly.
- `hla_mode: reads` calls directly from the aligned BAM and is faster, but is
  limited to HLA-A, HLA-B, and HLA-C.

## Run And Resume

Launch from the working directory containing the named run:

```bash
longpipe run 260410_BT-brain
```

Run only FIRE and its required upstream dependencies with:

```bash
longpipe run 260410_BT-brain --only-fire
```

This target skips cCRE aggregation, Sawfish, SVTopo, and HiFiHLA. If the
nucleosome-annotated or phased WGS inputs are missing from scratch, Snakemake
builds them because FIRE requires them. FIRE-only publication updates the
pileup, bigWigs, and haplotype analysis while preserving already published
PBWGS, SV, HLA, and cCRE outputs.

The command remains attached until the workflow succeeds or fails. All active
Snakemake state and intermediate files are written below:

```text
/hpc/temp/furlan_s/260410_BT-brain/
```

Longpipe snapshots the configuration and runsheet for every invocation.
Changes to the durable files cannot affect an active invocation. A later
`longpipe run` snapshots the current files and resumes the retained Snakemake
state.

To explicitly discard existing scratch state:

```bash
longpipe run 260410_BT-brain --restart
```

Concurrent invocations with the same run name are rejected. Scratch ownership
also prevents a run in another project directory from accidentally resuming or
restarting the same scratch path.

## Published Data

After a successful workflow, Longpipe stages and publishes:

```text
260410_BT-brain/
├── config.yaml
├── runsheet.csv
├── samples/
│   ├── Brain_1/
│   │   ├── pbWGS/
│   │   ├── fire/
│   │   │   ├── pileups/
│   │   │   ├── bigwigs/
│   │   │   ├── haplotypes/
│   │   │   │   ├── tracks/
│   │   │   │   └── differences/
│   │   │   └── ccre/
│   │   │       ├── encode/
│   │   │       └── brain/
│   │   ├── sawfish/
│   │   ├── svtopo/
│   │   └── hla/
│   └── Brain_2/
│       └── ...
├── logs/
├── benchmarks/
└── run-metadata/
    └── fire-eligibility.tsv
```

The `fire/` directory is present only for samples that enter the Fiber-seq
branch and pass the configured mean-depth threshold. Published FIRE pileups
and haplotype-difference tables include their Tabix `.tbi` indexes. The
haplotype difference table contains H1/H2 coverage, percent accessibility,
percentage-point difference, odds ratio, raw p-value, adjusted q-value, and a
significance flag; its companion summary records filtering and test counts.
Scratch paths retain the workflow's internal `results/{wgs,fire,sv,hla}`
layout so existing work can resume; Longpipe performs the per-sample
transformation only when publishing a successful run. A successful full
publication removes the superseded durable `results/` tree from older Longpipe
versions.

`run-metadata/` contains the latest status and versioned invocation records
with the exact configuration, runsheet, and selected Snakemake metadata used
for each launch.

On failure, Longpipe publishes logs, benchmarks, provenance, and failure status
but does not publish partial scientific results. Scratch is retained after both
success and failure.

Remove scratch explicitly after results have been checked if desired:

```bash
longpipe clean 260410_BT-brain
```

This never removes durable configuration or published results.

## Resource Tuning

Rules request total memory with `mem_mb`; the SLURM executor maps it to
`--mem`. Cluster-specific settings live in `profiles/slurm/config.yaml`. The
profile limits concurrent nested MiniWDL drivers, constrains large copy
operations, requests a GPU only for nucleosome annotation, and emits a SLURM
efficiency report. Up to ten sample-level MiniWDL workflows run at a time;
each `hifi_wgs` job consumes one `wdl_driver` resource unit. Native FIRE jobs
are independently scheduled and can resume at rule boundaries.

The cCRE rule performs one `bedtools map` pass for max, median, and sum. With
`scatac_sorted: true`, the cCRE BED must use the chromosome order in
`chrom_sizes`; Longpipe streams indexed FIRE pileups in that order.

## Development

Run the CLI tests with:

```bash
python -m unittest discover -s tests -v
```
