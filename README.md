# longpipe

`longpipe` runs a Snakemake workflow for PacBio HiFi and Fiber-seq samples. It
includes:

1. optional Fiber-seq m6A prediction or nucleosome annotation;
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
environment needs to be activated for standard m6A-tagged inputs.

Longpipe pins fibertools-rs `0.12.0` in both `envs/fibertools.yaml` and
`envs/fire-core.yaml`. Keeping prediction, nucleosome annotation, FIRE
calling, extraction, and pileup on one Fibertools release ensures that FIRE
reads the current `MA` nucleosome/MSP representation. The FIRE environment
also pins Samtools and HTSlib `1.23.1`; Mosdepth `0.3.7` requires HTSlib below
1.24, while fibertools-rs `0.12.0` requires the newer `libzlib` stack.

Kinetics inputs use `ft predict-m6a`. The portable default runs the pinned
Bioconda build on CPUs because Bioconda does not include libtorch or GPU
acceleration. This fallback requires no separate Fibertools installation, but
it is substantially slower. See the
[Fibertools installation guide](https://fiberseq.github.io/fibertools/install.html)
for the upstream explanation of the Bioconda and libtorch builds.

Users with a prebuilt GPU-capable executable can opt in per run:

```yaml
ft_predict_backend: external
ft_predict_executable: /absolute/path/to/ft
ft_predict_expected_version: 0.12.0
```

The executable may instead be a site wrapper that fully prepares its own
runtime. Longpipe resolves run-relative paths against the durable run
directory, but rejects PATH-only external commands. Before creating or
resuming scratch state, `longpipe run` executes `ft --version` and
`ft predict-m6a --help` with inherited libtorch variables removed. This catches
missing shared libraries and unusable wrappers before Snakemake submits jobs.
When `ft_predict_expected_version` is set, the detected version must match.

Existing configurations remain compatible: absent `ft_predict_backend`, a
plain `ft_predict_executable: ft` selects Conda, while any previously
configured custom executable selects the external backend.

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

Longpipe retries a MiniWDL invocation once when its log matches the known
call-cache coherence `FileNotFoundError` race. Other MiniWDL failures are not
retried. The deployment wrapper retains each attempt below the sample's
MiniWDL run directory as `longpipe-miniwdl-attempt-<n>.log`. This is a
defensive workflow measure; the cache race should also be reported and
corrected upstream in MiniWDL.

The deployment wrapper is shipped under `share/longpipe/deploy/miniwdl` in
the Longpipe installation. A Longpipe upgrade does not silently overwrite the
shared MiniWDL deployment. Update its outer `bin/miniwdl` wrapper deliberately;
the managed MiniWDL environment and cache contents do not need modification.

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
| `run` | PacBio run directory below `rawdir`; used for `raw` input when `hifi_bam` is blank |
| `cell` | SMRT Cell directory used for raw-BAM discovery |
| `barcode` | Barcode used for raw-BAM discovery |
| `hifi_bam` | Optional explicit raw BAM, or required recovery BAM for `aligned` input |
| `input_stage` | `raw` (default), `aligned`, or `pbwgs` |
| `pbwgs_outputs_json` | Existing PacBio WDL `outputs.json`; required only for `pbwgs` input |
| `sex` | `M`, `MALE`, `F`, `FEMALE`, or blank for `default_sex` |
| `hia5` | Explicit yes/no value controlling the Fiber-seq branch |
| `fiberseq_input` | Optional Fiber-seq input override: `auto`, `m6a`, or `kinetics`; blank defaults to automatic detection |
| `notes` | Optional free-text metadata retained with the runsheet and ignored by the workflow |

Labels must be unique and contain only letters, digits, `.`, `_`, and `-`.
Configuration and runsheet schemas are checked before jobs launch.

### Reanalysis Inputs

`input_stage` makes reuse of processed data explicit. It is never inferred
from a BAM filename:

- `raw` is the normal mode and the default for old runsheets. Longpipe locates
  a BAM from `run`, `cell`, and `barcode`, or uses an explicit `hifi_bam`.
- `aligned` is a recovery mode. Set `hifi_bam` to a coordinate-sorted aligned
  or haplotagged BAM. Longpipe validates it, then supplies its reads to the
  currently configured PacBio WDL, which realigns them and reruns the current
  analysis. Leave `pbwgs_outputs_json` blank.
- `pbwgs` reuses a completed PacBio WDL analysis. Set
  `pbwgs_outputs_json` and leave `hifi_bam`, `run`, `cell`, and `barcode`
  blank. Longpipe skips MiniWDL and uses the manifest's
  `humanwgs_singleton.merged_haplotagged_bam` output.

For example, rerun the current WDL from a prior aligned BAM with:

```csv
core_name,lab_name,run,cell,barcode,hifi_bam,input_stage,pbwgs_outputs_json,sex,hia5,fiberseq_input
LL1_S1,,,,,/path/to/LL1_S1.GRCh38.haplotagged.bam,aligned,,F,yes,auto
```

Reuse a completed PacBio WDL run with:

```csv
core_name,lab_name,run,cell,barcode,hifi_bam,input_stage,pbwgs_outputs_json,sex,hia5,fiberseq_input
LL1_S1,,,,,,pbwgs,/path/to/pacbiowdlR/outputs.json,F,yes,auto
```

The reused manifest may be MiniWDL's flat output object or an object containing
an `outputs` mapping. Its BAM must still exist and must be named
`<sample>.GRCh38.haplotagged.bam`, where `<sample>` is the runsheet label. For
a Fiber-seq sample, the manifest must also contain either the current
`humanwgs_singleton.stat_depth_mean` output or legacy
`humanwgs_singleton.stat_mean_depth` output.

A full `pbwgs` run copies and publishes the reused PacBio WDL results, then
runs Longpipe's FIRE, Sawfish, SVTopo, and HiFiHLA branches. A `pbwgs`
`--only-fire` run reads the existing manifest and BAM directly; it does not
copy the complete PacBio WDL result tree or rerun the other branches.

Reusing a version 2 PacBio WDL manifest does not update its variant calls to
version 3 methods. Use `pbwgs` when those existing calls are acceptable, such
as a FIRE-only analysis. Use `aligned` or, preferably, the original raw BAM
when the goal is a current complete WDL analysis.

Both recovery modes run `samtools quickcheck`, require coordinate sorting,
require an exact match to the configured reference sequence dictionary, and
confirm that sampled records include mapped primary reads. These checks establish
format and reference compatibility; they cannot prove that an aligned BAM
retains every original read or tag. Validation reports are published under
`run-metadata/reanalysis-input/`.

### Fiber-seq Input Modes

For `hia5: yes` samples, Longpipe inspects up to 5,000 primary reads before
selecting the preparation step. Blank or `auto` is recommended:

- Existing `A+a` calls in `MM/ML`, without nucleosome and MSP annotations,
  run `ft add-nucleosomes`. Longpipe recognizes both the current `MA` encoding
  (`nuc` and `msp` types) and legacy `ns`/`nl`/`as`/`al` tags.
- Existing m6A plus nucleosome and MSP annotations use the source BAM
  directly.
- No m6A plus PacBio `fi`, `ri`, `fp`, and `rp` kinetics runs
  `ft predict-m6a`. Prediction adds nucleosome and MSP annotations, so
  `ft add-nucleosomes` is not run afterward.
- A BAM containing neither sufficient m6A nor sufficient kinetics fails with
  an actionable error instead of entering the wrong branch.

Detection uses the same 5,000-read and 50% defaults as `ft validate`; these
are configurable with `fiberseq_detection_reads` and
`fiberseq_detection_min_fraction`. Explicit `m6a` and `kinetics` runsheet
values override branch selection but still validate that the required tags
are present. The observed counts, fractions, reason, and selected action are
published for each Fiber-seq sample under `run-metadata/fiberseq-input/`.

Prediction uses 30 threads and 128 GB total memory. The default Conda backend
requests no GPU and allows 24 hours for the slower CPU implementation. The
external backend requests one GPU and allows eight hours. Unlike
`--mem-per-cpu=21000MB -c 30`, both modes request 128 GB total rather than
approximately 630 GB. The rule records the selected backend and actual
Fibertools version in its log. `ft predict-m6a` drops the bulky kinetics tags
by default, matching the standard command; the source BAM is read in place and
remains unchanged.

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

This target skips cCRE aggregation, Sawfish, SVTopo, and HiFiHLA. For `raw` and
`aligned` inputs, Snakemake builds missing nucleosome-annotated or phased WGS
inputs because FIRE requires them. For `pbwgs`, it reuses the manifest's final
BAM instead. FIRE-only publication updates the pileup, bigWigs, and haplotype
analysis while preserving already published PBWGS, SV, HLA, and cCRE outputs.

Run only specific runsheet rows with `--rows`:

```bash
longpipe run 260410_BT-brain --rows 1,3,4:8
```

Rows are numbered from 1 after the CSV header, ranges are inclusive, and the
example selects rows 1, 3, 4, 5, 6, 7, and 8. Commas combine selections;
duplicate selections are harmless. Invalid, reversed, zero, and out-of-range
rows are rejected before the workflow starts. The durable `runsheet.csv`
remains unchanged, while the invocation record contains the filtered runsheet
actually used. Successful subset runs replace results only for the selected
samples and preserve published results and metadata for all unselected
samples. `--rows` can also be combined with `--only-fire`.

The command remains attached until the workflow succeeds or fails. All active
Snakemake state and intermediate files are written below:

```text
/hpc/temp/furlan_s/260410_BT-brain/
```

Longpipe snapshots the configuration and runsheet for every invocation.
Changes to the durable files cannot affect an active invocation. A later
`longpipe run` snapshots the current files and resumes the retained Snakemake
state.

Preview the exact workflow that Longpipe would resume without submitting jobs
or publishing results:

```bash
longpipe run 260410_BT-brain --dry-run
```

The dry-run uses the normal configuration snapshot, run lock, scratch state,
SLURM profile, target selection, and optional `--rows`/`--only-fire` filters.
Its invocation and status are recorded under `run-metadata/`, but scientific
results and diagnostics are not published. `--dry-run` cannot be combined with
`--restart`, because a preview must not discard retained scratch state.

To explicitly discard existing scratch state:

```bash
longpipe run 260410_BT-brain --restart
```

Concurrent invocations with the same run name are rejected. Scratch ownership
also prevents a run in another project directory from accidentally resuming or
restarting the same scratch path.

## Published Data

Each sample follows its own Snakemake trajectory. As soon as all outputs for
one sample are complete, Longpipe atomically publishes that sample without
waiting for other samples in the run:

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
    ├── fire-eligibility.tsv
    ├── sample-publication/
    │   └── <sample>.json
    └── fire-publication/
        └── <sample>.json
```

The `fire/` directory is present only for samples that enter the Fiber-seq
branch and pass the configured mean-depth threshold. Published FIRE pileups
and haplotype-difference tables include their Tabix `.tbi` indexes. The
haplotype difference table contains H1/H2 coverage, percent accessibility,
percentage-point difference, odds ratio, raw p-value, adjusted q-value, and a
significance flag; its companion summary records filtering and test counts.
Scratch paths retain the workflow's internal `results/{wgs,fire,sv,hla}`
layout so existing work can resume. Longpipe performs the layout transformation
independently when each sample reaches its publication rule. A successful full
run removes the superseded durable `results/` tree from older Longpipe versions.

`run-metadata/` contains the latest status and versioned invocation records
with the exact configuration, runsheet, and selected Snakemake metadata used
for each launch. `run-metadata/fiberseq-input/` records the BAM tag evidence
and preparation action selected for every Fiber-seq sample.
`run-metadata/reanalysis-input/` records validation and provenance for every
`aligned` or `pbwgs` input.
`run-metadata/sample-publication/` records which full sample trees completed
publication. FIRE-only runs record their independent sample-level updates in
`run-metadata/fire-publication/`. The combined `fire-eligibility.tsv` is a
run-level summary and is finalized after every Fiber-seq candidate has an
eligibility decision; it does not gate an individual sample's analysis.

On failure, Longpipe publishes logs, benchmarks, provenance, and failure status.
Any sample that already reached its atomic publication rule remains available
as a complete sample result; incomplete samples are not published. Scratch is
retained after both success and failure.

Remove scratch explicitly after results have been checked if desired:

```bash
longpipe clean 260410_BT-brain
```

This never removes durable configuration or published results.

## Resource Tuning

Rules request total memory with `mem_mb`; the SLURM executor maps it to
`--mem`. Cluster-specific settings live in `profiles/slurm/config.yaml`. The
profile limits concurrent nested MiniWDL drivers, constrains large copy
operations, requests a GPU only for m6A prediction, and emits a SLURM
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
