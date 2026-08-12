import glob
import json
import math
import os
import re
from pathlib import Path

import pandas as pd
from snakemake.utils import min_version, validate


min_version("8.6")

envvars:
    "LONGPIPE_RUN_DIR",
    "LONGPIPE_RUN_NAME",
    "LONGPIPE_RUNSHEET",
    "LONGPIPE_SCRATCH_DIR"

WORKFLOW_ROOT = Path(str(workflow.source_path("Snakefile"))).resolve().parent
CONFIG_SCHEMA = str(workflow.source_path("schemas/config.schema.yaml"))
RUNSHEET_SCHEMA = str(workflow.source_path("schemas/runsheet.schema.yaml"))

validate(config, CONFIG_SCHEMA)

RUN_DIRECTORY = Path(os.environ["LONGPIPE_RUN_DIR"]).resolve()


def configured_path(value):
    path = Path(os.path.expandvars(os.path.expanduser(str(value))))
    if not path.is_absolute():
        path = RUN_DIRECTORY / path
    return str(path.resolve())


def configured_executable(value):
    executable = os.path.expandvars(os.path.expanduser(str(value)))
    if os.sep not in executable:
        return executable
    return configured_path(executable)


NAME = os.environ["LONGPIPE_RUN_NAME"]
OUTDIR = str(Path(os.environ["LONGPIPE_SCRATCH_DIR"]).resolve())
RAWDIR = configured_path(config["rawdir"])
RUNSHEET = str(Path(os.environ["LONGPIPE_RUNSHEET"]).resolve())

HIA5_GATES = config["hia5_gates_fiberseq"]
HLA_MODE = config["hla_mode"]
SCATAC_SORTED = config["scatac_sorted"]
FIRE_VERSION = config["fire_version"]
FIRE_MIN_MEAN_DEPTH = float(config.get("fire_min_mean_depth", 8.0))
FIRE_UPSTREAM_COMMIT = config["fire_upstream_commit"]
FIRE_MIN_COVERAGE = int(config.get("fire_min_coverage", 4))
FIRE_COVERAGE_WITHIN_N_SD = float(
    config.get("fire_coverage_within_n_sd", 5.0)
)
FIRE_FILTER_FLAG = int(config.get("fire_filter_flag", 260))
FIRE_MIN_MSP = int(config.get("fire_min_msp", 10))
FIRE_MIN_AVERAGE_MSP_SIZE = int(
    config.get("fire_min_average_msp_size", 10)
)
FIRE_MIN_CONTIG_LENGTH = int(config.get("fire_min_contig_length", 0))
FIRE_KEEP_CHROMOSOMES = config.get("fire_keep_chromosomes", ".*")
FIRE_HAPLOTYPE_MIN_COVERAGE = int(
    config.get("fire_haplotype_min_coverage", 10)
)
FIRE_HAPLOTYPE_MAX_FDR = float(config.get("fire_haplotype_max_fdr", 0.05))
FIRE_HAPLOTYPE_QVALUE = float(config.get("fire_haplotype_qvalue", 0.05))
FIRE_MEAN_DEPTH_OUTPUT_KEYS = (
    "humanwgs_singleton.stat_depth_mean",
    "humanwgs_singleton.stat_mean_depth",
)

FASTA = configured_path(config["ref_fasta"])
FIRE_REF = configured_path(config["fire_ref"])
HLA_REF = configured_path(config["hla_ref_fasta"])
CHROM_SIZES = configured_path(config["chrom_sizes"])
CCRE_PEAKS = {
    name: configured_path(path)
    for name, path in config["ccre_peaks"].items()
}
CCRE_PEAK_SETS = list(CCRE_PEAKS)
REF_MAP = configured_path(config["ref_map_file"])
TERTIARY_MAP = configured_path(config["tertiary_map_file"])
WDL = configured_path(config["wdl_singleton"])
FIRE_EXCLUDE_REGIONS = [
    configured_path(path) for path in config["fire_exclude_regions"]
]
CNV_EXCLUDED = configured_path(config["cnv_excluded_regions"])
EXPECTED_CN = {
    "MALE": configured_path(config["expected_cn_xy"]),
    "FEMALE": configured_path(config["expected_cn_xx"]),
}

COPY_PBWDL = configured_path(config["copy_pbwdl_tool"])
SHRINK_FIRE_BED = configured_path(config["shrink_fire_bed_tool"])
MINIWDL = configured_executable(config["miniwdl_executable"])
FT_PREDICT_CONFIGURED = config.get("ft_predict_executable", "ft")
FT_PREDICT_BACKEND = config.get("ft_predict_backend")
if FT_PREDICT_BACKEND is None:
    FT_PREDICT_BACKEND = (
        "external" if FT_PREDICT_CONFIGURED != "ft" else "conda"
    )
FT_PREDICT = (
    configured_executable(FT_PREDICT_CONFIGURED)
    if FT_PREDICT_BACKEND == "external"
    else "ft"
)
FIBERSEQ_DETECTION_READS = int(config.get("fiberseq_detection_reads", 5000))
FIBERSEQ_DETECTION_MIN_FRACTION = float(
    config.get("fiberseq_detection_min_fraction", 0.5)
)

FIRE_PILEUP = (
    f"results/fire/pileups/"
    f"{{sample}}-fire-{FIRE_VERSION}-pileup.bed.gz"
)
FIRE_PILEUP_INDEX = f"{FIRE_PILEUP}.tbi"
FIRE_BIGWIG = (
    f"results/fire/bigwigs/"
    f"{{sample}}-fire-{FIRE_VERSION}-pileup_shrunk.bw"
)
FIRE_H1_BIGWIG = (
    f"results/fire/haplotypes/tracks/"
    f"{{sample}}-fire-{FIRE_VERSION}-H1-percent-accessible.bw"
)
FIRE_H2_BIGWIG = (
    f"results/fire/haplotypes/tracks/"
    f"{{sample}}-fire-{FIRE_VERSION}-H2-percent-accessible.bw"
)
FIRE_HAPLOTYPE_DIFFERENCE_BIGWIG = (
    f"results/fire/haplotypes/tracks/"
    f"{{sample}}-fire-{FIRE_VERSION}-H1-minus-H2-percent-accessible.bw"
)
FIRE_HAPLOTYPE_DIFFERENCES = (
    f"results/fire/haplotypes/differences/"
    f"{{sample}}-fire-{FIRE_VERSION}-haplotype-differences.bed.gz"
)
FIRE_HAPLOTYPE_DIFFERENCES_INDEX = f"{FIRE_HAPLOTYPE_DIFFERENCES}.tbi"
FIRE_HAPLOTYPE_SUMMARY = (
    f"results/fire/haplotypes/differences/"
    f"{{sample}}-fire-{FIRE_VERSION}-haplotype-differences-summary.tsv"
)
FIRE_ELIGIBILITY = "intermediate/fire/eligibility/{sample}.json"
FIRE_ELIGIBILITY_ROW = "intermediate/fire/eligibility/{sample}.tsv"
FIRE_ELIGIBILITY_REPORT = "results/fire/eligibility.tsv"
PUBLISHED_FIRE_ELIGIBILITY_REPORT = str(
    RUN_DIRECTORY / "run-metadata" / "fire-eligibility.tsv"
)
SAMPLE_PUBLICATION_STATUS = str(
    RUN_DIRECTORY
    / "run-metadata"
    / "sample-publication"
    / "{sample}.json"
)
FIRE_PUBLICATION_STATUS = str(
    RUN_DIRECTORY
    / "run-metadata"
    / "fire-publication"
    / "{sample}.json"
)
FIBERSEQ_DETECTION_REPORT = (
    "intermediate/fiberseq/input-detection/{sample}.json"
)
REANALYSIS_VALIDATION_REPORT = (
    "intermediate/reanalysis/input-validation/{sample}.json"
)
PBWGS_OUTPUTS = "pbWGS/{sample}/outputs.json"
HAPLOTAGGED_BAM = (
    "results/wgs/{sample}/{sample}.GRCh38.haplotagged.bam"
)
PUBLISHED_WGS_MANIFEST = "results/wgs/{sample}/outputs.json"
HAPLOTAGGED_BAM_INDEX = f"{HAPLOTAGGED_BAM}.bai"
FIRE_INPUT_BAM = "intermediate/fire/input/{sample}.bam"
FIRE_INPUT_BAM_INDEX = f"{FIRE_INPUT_BAM}.bai"
HLA_JSON = "results/hla/{sample}_hifihla_report.json"
HLA_REPORT = "results/hla/{sample}_hifihla_report.tsv"
HLA_SUMMARY = "results/hla/{sample}_hifihla_summary.tsv"


required_columns = {
    "core_name",
    "lab_name",
    "run",
    "cell",
    "barcode",
    "sex",
    "hia5",
}
sample_table = pd.read_csv(RUNSHEET, dtype=str, keep_default_na=False)
missing_columns = sorted(required_columns - set(sample_table.columns))
if missing_columns:
    raise ValueError(
        f"{RUNSHEET}: missing required columns: {', '.join(missing_columns)}"
    )

validate(sample_table, RUNSHEET_SCHEMA)


def sample_label(row):
    lab_name = row["lab_name"].strip()
    core_name = row["core_name"].strip()
    if not (lab_name or core_name):
        raise ValueError("each runsheet row needs core_name or lab_name")
    return lab_name or core_name


sample_table["label"] = sample_table.apply(sample_label, axis=1)
invalid_labels = [
    label
    for label in sample_table["label"]
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", label)
]
if invalid_labels:
    raise ValueError(
        "sample labels may contain only letters, digits, '.', '_', and '-': "
        + ", ".join(invalid_labels)
    )

sample_table = sample_table.set_index("label", drop=False)
duplicate_labels = sorted(
    set(sample_table.index[sample_table.index.duplicated()].tolist())
)
if duplicate_labels:
    raise ValueError(f"duplicate sample labels: {', '.join(duplicate_labels)}")
if sample_table.empty:
    raise ValueError(f"{RUNSHEET}: no samples found")

SAMPLES = sample_table.index.tolist()


def parse_boolean(value, field, sample):
    normalized = str(value).strip().lower()
    if normalized in {"y", "yes", "true", "t", "1"}:
        return True
    if normalized in {"n", "no", "false", "f", "0"}:
        return False
    raise ValueError(f"{sample}: invalid {field} value {value!r}")


def is_hia5(sample):
    return parse_boolean(sample_table.loc[sample, "hia5"], "hia5", sample)


def uses_fiberseq(sample):
    return not HIA5_GATES or is_hia5(sample)


def input_stage_of(sample):
    value = sample_table.loc[sample].get("input_stage", "").strip().lower()
    if not value:
        return "raw"
    if value not in {"raw", "aligned", "pbwgs"}:
        raise ValueError(
            f"{sample}: input_stage must be raw, aligned, pbwgs, or blank; "
            f"got {value!r}"
        )
    return value


def fiberseq_input_of(sample):
    value = sample_table.loc[sample].get("fiberseq_input", "").strip().lower()
    if not value:
        return "auto"
    if value not in {"auto", "m6a", "kinetics"}:
        raise ValueError(
            f"{sample}: fiberseq_input must be auto, m6a, kinetics, or blank; "
            f"got {value!r}"
        )
    return value


def sex_of(sample):
    value = sample_table.loc[sample, "sex"].strip().upper()
    if not value:
        value = config["default_sex"]
    aliases = {"M": "MALE", "MALE": "MALE", "F": "FEMALE", "FEMALE": "FEMALE"}
    if value not in aliases:
        raise ValueError(
            f"{sample}: sex must be M, MALE, F, FEMALE, or blank; got {value!r}"
        )
    return aliases[value]


def hifi_bam(sample):
    row = sample_table.loc[sample]
    explicit_bam = row.get("hifi_bam", "").strip()
    if explicit_bam:
        bam_path = configured_path(explicit_bam)
        if not os.path.isfile(bam_path):
            raise ValueError(f"{sample}: hifi_bam does not exist: {bam_path}")
        return bam_path

    missing_fields = [
        field for field in ("run", "cell", "barcode") if not row[field].strip()
    ]
    if missing_fields:
        raise ValueError(
            f"{sample}: set hifi_bam or provide {', '.join(missing_fields)}"
        )

    barcode = re.sub(r"^bc-", "bc", row["barcode"].strip(), flags=re.IGNORECASE)
    pattern = os.path.join(
        RAWDIR,
        row["run"].strip(),
        row["cell"].strip(),
        "hifi_reads",
        f"*hifi_reads.{barcode}.bam",
    )
    matches = sorted(glob.glob(pattern))
    if len(matches) != 1:
        raise ValueError(
            f"{sample}: expected exactly one HiFi BAM matching {pattern}; "
            f"found {matches}"
        )
    return str(Path(matches[0]).resolve())


_existing_pbwgs_cache = {}


def existing_pbwgs(sample):
    if sample in _existing_pbwgs_cache:
        return _existing_pbwgs_cache[sample]

    configured = sample_table.loc[sample].get("pbwgs_outputs_json", "").strip()
    if not configured:
        raise ValueError(
            f"{sample}: input_stage=pbwgs requires pbwgs_outputs_json"
        )
    manifest = Path(configured_path(configured))
    if not manifest.is_file():
        raise ValueError(
            f"{sample}: pbwgs_outputs_json does not exist: {manifest}"
        )
    try:
        with manifest.open(encoding="utf-8") as handle:
            payload = json.load(handle)
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(
            f"{sample}: could not read pbwgs_outputs_json {manifest}: {error}"
        ) from error
    if not isinstance(payload, dict):
        raise ValueError(
            f"{sample}: pbwgs_outputs_json is not a JSON object: {manifest}"
        )
    outputs = payload.get("outputs", payload)
    if not isinstance(outputs, dict):
        raise ValueError(
            f"{sample}: pbwgs_outputs_json outputs are not a JSON object: "
            f"{manifest}"
        )
    bam_value = outputs.get("humanwgs_singleton.merged_haplotagged_bam")
    if not isinstance(bam_value, str) or not bam_value:
        raise ValueError(
            f"{sample}: {manifest} does not contain "
            "humanwgs_singleton.merged_haplotagged_bam"
        )
    bam = Path(os.path.expandvars(os.path.expanduser(bam_value)))
    if not bam.is_absolute():
        bam = manifest.parent / bam
    bam = bam.resolve()
    if not bam.is_file():
        raise ValueError(
            f"{sample}: reused PacBio WDL BAM does not exist: {bam}"
        )
    expected_name = f"{sample}.GRCh38.haplotagged.bam"
    if bam.name != expected_name:
        raise ValueError(
            f"{sample}: reused PacBio WDL BAM must be named {expected_name}; "
            f"found {bam.name}. Use the original PacBio sample label."
        )

    result = {
        "manifest": str(manifest.resolve()),
        "bam": str(bam),
        "workflow_version": outputs.get("humanwgs_singleton.workflow_version"),
    }
    _existing_pbwgs_cache[sample] = result
    return result


for _sample in SAMPLES:
    _stage = input_stage_of(_sample)
    _row = sample_table.loc[_sample]
    _explicit_bam = _row.get("hifi_bam", "").strip()
    _outputs_json = _row.get("pbwgs_outputs_json", "").strip()
    if _stage == "aligned" and not _explicit_bam:
        raise ValueError(
            f"{_sample}: input_stage=aligned requires an explicit hifi_bam"
        )
    if _stage == "pbwgs":
        if _explicit_bam:
            raise ValueError(
                f"{_sample}: input_stage=pbwgs cannot also set hifi_bam"
            )
        existing_pbwgs(_sample)
    elif _outputs_json:
        raise ValueError(
            f"{_sample}: pbwgs_outputs_json requires input_stage=pbwgs"
        )


def source_bam(sample):
    if input_stage_of(sample) == "pbwgs":
        return existing_pbwgs(sample)["bam"]
    return hifi_bam(sample)


def pbwgs_outputs_for_sample(sample):
    if input_stage_of(sample) == "pbwgs":
        return existing_pbwgs(sample)["manifest"]
    return PBWGS_OUTPUTS.format(sample=sample)


def reanalysis_validation_input(wildcards):
    if input_stage_of(wildcards.sample) == "raw":
        return []
    return REANALYSIS_VALIDATION_REPORT.format(sample=wildcards.sample)


def fiberseq_detection_checkpoint(wildcards):
    return checkpoints.detect_fiberseq_input.get(sample=wildcards.sample)


def wgs_hifi(wildcards):
    if input_stage_of(wildcards.sample) == "pbwgs":
        raise ValueError(
            f"{wildcards.sample}: reused pbWGS outputs must not enter hifi_wgs"
        )
    if not uses_fiberseq(wildcards.sample):
        return hifi_bam(wildcards.sample)

    checkpoint_output = fiberseq_detection_checkpoint(wildcards)
    with open(checkpoint_output.output.report, encoding="utf-8") as handle:
        action = json.load(handle)["selected_action"]
    if action == "use_existing":
        return hifi_bam(wildcards.sample)
    if action == "add_nucleosomes":
        return f"fibertools/{wildcards.sample}.nucleosomes.bam"
    if action == "predict_m6a":
        return f"fibertools/{wildcards.sample}.m6a_predict.bam"
    raise ValueError(
        f"{wildcards.sample}: unrecognized Fiber-seq input action {action!r}"
    )


def prepared_external_fiberseq_bam(wildcards):
    checkpoint_output = fiberseq_detection_checkpoint(wildcards)
    with open(checkpoint_output.output.report, encoding="utf-8") as handle:
        action = json.load(handle)["selected_action"]
    if action == "use_existing":
        return source_bam(wildcards.sample)
    if action == "add_nucleosomes":
        return f"fibertools/{wildcards.sample}.nucleosomes.bam"
    if action == "predict_m6a":
        return f"fibertools/{wildcards.sample}.m6a_predict.bam"
    raise ValueError(
        f"{wildcards.sample}: unrecognized Fiber-seq input action {action!r}"
    )


def fire_source_bam(wildcards):
    if input_stage_of(wildcards.sample) == "pbwgs":
        return prepared_external_fiberseq_bam(wildcards)
    return HAPLOTAGGED_BAM.format(sample=wildcards.sample)


def svtopo_prefix(sample):
    return sample.replace("_", "-")


svtopo_prefixes = [svtopo_prefix(sample) for sample in SAMPLES]
if len(svtopo_prefixes) != len(set(svtopo_prefixes)):
    raise ValueError(
        "sample labels collide after replacing underscores for SVTopo prefixes"
    )

FIRE_CANDIDATES = [sample for sample in SAMPLES if uses_fiberseq(sample)]


def load_fire_chromosomes():
    if not FIRE_CANDIDATES:
        return []
    reference_index = Path(f"{FIRE_REF}.fai")
    if not reference_index.is_file():
        raise ValueError(f"FIRE reference index does not exist: {reference_index}")
    try:
        keep_pattern = re.compile(FIRE_KEEP_CHROMOSOMES)
    except re.error as error:
        raise ValueError(
            f"invalid fire_keep_chromosomes pattern: {FIRE_KEEP_CHROMOSOMES!r}"
        ) from error

    chromosomes = []
    with reference_index.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(
                    f"{reference_index}:{line_number}: malformed FASTA index"
                )
            chromosome = fields[0]
            try:
                length = int(fields[1])
            except ValueError as error:
                raise ValueError(
                    f"{reference_index}:{line_number}: invalid contig length"
                ) from error
            if "/" in chromosome:
                raise ValueError(
                    f"{reference_index}:{line_number}: contig contains '/': "
                    f"{chromosome}"
                )
            if length < FIRE_MIN_CONTIG_LENGTH:
                continue
            if "chrUn_" in chromosome or "_random" in chromosome:
                continue
            if keep_pattern.fullmatch(chromosome):
                chromosomes.append(chromosome)

    chromosomes.sort()
    if not chromosomes:
        raise ValueError(
            "no FIRE chromosomes remain after applying contig filters"
        )
    return chromosomes


FIRE_CHROMS = load_fire_chromosomes()
FIRE_CHROM_PATTERN = (
    "(?:" + "|".join(map(re.escape, FIRE_CHROMS)) + ")"
    if FIRE_CHROMS
    else r"(?!)"
)


def fire_eligibility_checkpoint(wildcards):
    return checkpoints.fire_eligibility.get(sample=wildcards.sample)


def fire_eligible(wildcards):
    checkpoint_output = fire_eligibility_checkpoint(wildcards)
    with open(checkpoint_output.output.manifest, encoding="utf-8") as handle:
        decision = json.load(handle)
    if decision.get("sample") != wildcards.sample:
        raise ValueError(
            f"{checkpoint_output.output.manifest}: FIRE eligibility sample does "
            f"not match {wildcards.sample}"
        )
    eligible = decision.get("eligible")
    if not isinstance(eligible, bool):
        raise ValueError(
            f"{checkpoint_output.output.manifest}: FIRE eligibility decision "
            "must be boolean"
        )
    return eligible


def eligible_fire_bam(wildcards):
    if not fire_eligible(wildcards):
        raise ValueError(
            f"{wildcards.sample}: FIRE requires mean depth >= "
            f"{FIRE_MIN_MEAN_DEPTH:g}"
        )
    return FIRE_INPUT_BAM.format(sample=wildcards.sample)


def fire_core_outputs(sample):
    return (
        expand(FIRE_PILEUP, sample=sample)
        + expand(FIRE_PILEUP_INDEX, sample=sample)
        + expand(FIRE_BIGWIG, sample=sample)
        + expand(FIRE_H1_BIGWIG, sample=sample)
        + expand(FIRE_H2_BIGWIG, sample=sample)
        + expand(FIRE_HAPLOTYPE_DIFFERENCE_BIGWIG, sample=sample)
        + expand(FIRE_HAPLOTYPE_DIFFERENCES, sample=sample)
        + expand(FIRE_HAPLOTYPE_DIFFERENCES_INDEX, sample=sample)
        + expand(FIRE_HAPLOTYPE_SUMMARY, sample=sample)
    )


def fire_core_targets(wildcards):
    checkpoint_output = fire_eligibility_checkpoint(wildcards)
    targets = [
        checkpoint_output.output.manifest,
        checkpoint_output.output.report,
    ]
    if fire_eligible(wildcards):
        targets.extend(fire_core_outputs(wildcards.sample))
    return targets


def fire_targets(wildcards):
    targets = fire_core_targets(wildcards)
    if fire_eligible(wildcards):
        targets.extend(
            expand(
                "results/fire/ccre/{peak_set}/{sample}_scATAC_max_scores.bed",
                sample=wildcards.sample,
                peak_set=CCRE_PEAK_SETS,
            )
        )
        targets.extend(
            expand(
                "results/fire/ccre/{peak_set}/{sample}_scATAC_median_scores.bed",
                sample=wildcards.sample,
                peak_set=CCRE_PEAK_SETS,
            )
        )
        targets.extend(
            expand(
                "results/fire/ccre/{peak_set}/{sample}_scATAC_sum_scores.bed",
                sample=wildcards.sample,
                peak_set=CCRE_PEAK_SETS,
            )
        )
    return targets


def sample_fire_publication_inputs(wildcards):
    if wildcards.sample not in FIRE_CANDIDATES:
        return []
    return fire_targets(wildcards)


def fire_eligibility_rows(wildcards):
    return [
        checkpoints.fire_eligibility.get(sample=sample).output.report
        for sample in FIRE_CANDIDATES
    ]


workdir: OUTDIR

shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; ")

wildcard_constraints:
    sample="[A-Za-z0-9][A-Za-z0-9._-]*",
    peak_set="[A-Za-z0-9][A-Za-z0-9._-]*",
    fire_chrom=FIRE_CHROM_PATTERN

localrules: all, fire_only, fire_eligibility_report, wdl_input, hifi_wgs, \
    fire_eligibility


rule fire_only:
    input:
        samples=expand(FIRE_PUBLICATION_STATUS, sample=FIRE_CANDIDATES),
        eligibility=PUBLISHED_FIRE_ELIGIBILITY_REPORT


rule all:
    default_target: True
    input:
        samples=expand(SAMPLE_PUBLICATION_STATUS, sample=SAMPLES),
        eligibility=PUBLISHED_FIRE_ELIGIBILITY_REPORT


rule publish_sample_results:
    input:
        wgs_bam=HAPLOTAGGED_BAM,
        wgs_manifest=PUBLISHED_WGS_MANIFEST,
        fire=sample_fire_publication_inputs,
        sawfish_vcf="results/sv/{sample}.sawfish.vcf.gz",
        sawfish_index="results/sv/{sample}.sawfish.vcf.gz.tbi",
        sawfish_reads="results/sv/{sample}.supporting_reads.json.gz",
        svtopo="results/sv/{sample}.svtopo",
        hla_json=HLA_JSON,
        hla_report=HLA_REPORT,
        hla_summary=HLA_SUMMARY
    output:
        status=SAMPLE_PUBLICATION_STATUS
    threads: 1
    resources:
        mem_mb=4000,
        runtime=240,
        io_heavy=1
    params:
        mode="full",
        scratch_directory=OUTDIR,
        run_directory=str(RUN_DIRECTORY),
        fire_version=FIRE_VERSION,
        peak_sets=CCRE_PEAK_SETS,
        fire_decision=lambda wildcards: (
            FIRE_ELIGIBILITY.format(sample=wildcards.sample)
            if wildcards.sample in FIRE_CANDIDATES
            else ""
        ),
        invocation_id=os.environ.get("LONGPIPE_INVOCATION_ID", "")
    log:
        "logs/publish_sample_results/{sample}.log"
    benchmark:
        "benchmarks/publish_sample_results/{sample}.tsv"
    script:
        "workflow/scripts/publish_sample.py"


rule publish_fire_sample:
    input:
        fire_core_targets
    output:
        status=FIRE_PUBLICATION_STATUS
    threads: 1
    resources:
        mem_mb=4000,
        runtime=120,
        io_heavy=1
    params:
        mode="fire_only",
        scratch_directory=OUTDIR,
        run_directory=str(RUN_DIRECTORY),
        fire_version=FIRE_VERSION,
        peak_sets=[],
        fire_decision=lambda wildcards: FIRE_ELIGIBILITY.format(
            sample=wildcards.sample
        ),
        invocation_id=os.environ.get("LONGPIPE_INVOCATION_ID", "")
    log:
        "logs/publish_fire_sample/{sample}.log"
    benchmark:
        "benchmarks/publish_fire_sample/{sample}.tsv"
    script:
        "workflow/scripts/publish_sample.py"


rule validate_reanalysis_input:
    input:
        bam=lambda wildcards: source_bam(wildcards.sample),
        reference_index=f"{FASTA}.fai"
    output:
        report=REANALYSIS_VALIDATION_REPORT
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30
    params:
        stage=lambda wildcards: input_stage_of(wildcards.sample),
        outputs_json=lambda wildcards: (
            existing_pbwgs(wildcards.sample)["manifest"]
            if input_stage_of(wildcards.sample) == "pbwgs"
            else ""
        ),
        script=str(
            workflow.source_path(
                "workflow/scripts/validate_reanalysis_input.py"
            )
        )
    log:
        "logs/validate_reanalysis_input/{sample}.log"
    benchmark:
        "benchmarks/validate_reanalysis_input/{sample}.tsv"
    conda:
        "envs/fibertools.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.report:q})" "$(dirname {log:q})"
        outputs_json={params.outputs_json:q}
        run_validator() {{
            python {params.script:q} \
                --input {input.bam:q} \
                --output {output.report:q} \
                --reference-index {input.reference_index:q} \
                --sample {wildcards.sample:q} \
                --stage {params.stage:q} \
                "$@"
        }}
        if [[ -n "$outputs_json" ]]; then
            run_validator --outputs-json "$outputs_json"
        else
            run_validator
        fi > {log:q} 2>&1
        """


checkpoint detect_fiberseq_input:
    input:
        bam=lambda wildcards: source_bam(wildcards.sample),
        validation=reanalysis_validation_input
    output:
        report=FIBERSEQ_DETECTION_REPORT
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30
    params:
        mode=lambda wildcards: fiberseq_input_of(wildcards.sample),
        reads=FIBERSEQ_DETECTION_READS,
        min_fraction=FIBERSEQ_DETECTION_MIN_FRACTION,
        script=str(
            workflow.source_path(
                "workflow/scripts/detect_fiberseq_input.py"
            )
        )
    log:
        "logs/detect_fiberseq_input/{sample}.log"
    benchmark:
        "benchmarks/detect_fiberseq_input/{sample}.tsv"
    conda:
        "envs/fibertools.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.report:q})" "$(dirname {log:q})"
        python {params.script:q} \
            --input {input.bam:q} \
            --output {output.report:q} \
            --sample {wildcards.sample:q} \
            --mode {params.mode:q} \
            --reads {params.reads} \
            --min-fraction {params.min_fraction} \
            > {log:q} 2>&1
        """


rule add_nucleosomes:
    input:
        bam=lambda wildcards: source_bam(wildcards.sample)
    output:
        "fibertools/{sample}.nucleosomes.bam"
    threads: 16
    resources:
        mem_mb=64000,
        runtime=240
    log:
        "logs/add_nucleosomes/{sample}.log"
    benchmark:
        "benchmarks/add_nucleosomes/{sample}.tsv"
    conda:
        "envs/fibertools.yaml"
    shell:
        """
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        ft add-nucleosomes -t {threads} {input.bam:q} > {output:q} 2> {log:q}
        test -s {output:q}
        """


rule predict_m6a:
    input:
        bam=lambda wildcards: source_bam(wildcards.sample)
    output:
        "fibertools/{sample}.m6a_predict.bam"
    threads: 30
    resources:
        mem_mb=128000,
        runtime=480 if FT_PREDICT_BACKEND == "external" else 1440,
        gpu=1 if FT_PREDICT_BACKEND == "external" else 0
    params:
        executable=FT_PREDICT,
        backend=FT_PREDICT_BACKEND
    log:
        "logs/predict_m6a/{sample}.log"
    benchmark:
        "benchmarks/predict_m6a/{sample}.tsv"
    conda:
        "envs/fibertools.yaml"
    shell:
        """
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        printf 'backend: %s\n' {params.backend:q} > {log:q}
        {params.executable:q} --version >> {log:q} 2>&1
        {params.executable:q} predict-m6a -t {threads} \
            {input.bam:q} {output:q} >> {log:q} 2>&1
        test -s {output:q}
        """


rule wdl_input:
    input:
        bam=wgs_hifi,
        ref_map=REF_MAP,
        tertiary_map=TERTIARY_MAP,
        wdl=WDL
    output:
        "pbWGS/{sample}/input.json"
    run:
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "humanwgs_singleton.sample_id": wildcards.sample,
            "humanwgs_singleton.sex": sex_of(wildcards.sample),
            "humanwgs_singleton.hifi_reads": [str(Path(input.bam).resolve())],
            "humanwgs_singleton.ref_map_file": str(Path(input.ref_map).resolve()),
            "humanwgs_singleton.tertiary_map_file": str(
                Path(input.tertiary_map).resolve()
            ),
            "humanwgs_singleton.backend": "HPC",
            "humanwgs_singleton.preemptible": config["wdl_preemptible"],
        }
        with open(output[0], "w", encoding="utf-8") as output_handle:
            json.dump(payload, output_handle, indent=2)
            output_handle.write("\n")


rule hifi_wgs:
    input:
        bam=wgs_hifi,
        config="pbWGS/{sample}/input.json",
        wdl=WDL,
        validation=reanalysis_validation_input
    output:
        PBWGS_OUTPUTS
    threads: 1
    resources:
        mem_mb=4000,
        runtime=10080,
        wdl_driver=1
    params:
        executable=MINIWDL,
        run_dir="pbWGS/{sample}/runs"
    log:
        "logs/hifi_wgs/{sample}.log"
    benchmark:
        "benchmarks/hifi_wgs/{sample}.tsv"
    shell:
        """
        mkdir -p {params.run_dir:q} "$(dirname {output:q})" "$(dirname {log:q})"
        {params.executable:q} run {input.wdl:q} \
            -i {input.config:q} \
            --dir {params.run_dir:q} \
            -o {output:q} \
            --verbose > {log:q} 2>&1
        test -s {output:q}
        """


checkpoint fire_eligibility:
    input:
        manifest=lambda wildcards: pbwgs_outputs_for_sample(
            wildcards.sample
        ),
        validation=reanalysis_validation_input
    output:
        manifest=FIRE_ELIGIBILITY,
        report=FIRE_ELIGIBILITY_ROW
    params:
        minimum_mean_depth=FIRE_MIN_MEAN_DEPTH
    run:
        Path(output.manifest).parent.mkdir(parents=True, exist_ok=True)
        Path(output.report).parent.mkdir(parents=True, exist_ok=True)
        with open(input.manifest, encoding="utf-8") as handle:
            payload = json.load(handle)
        outputs = payload.get("outputs", payload)
        mean_depth = None
        for output_key in FIRE_MEAN_DEPTH_OUTPUT_KEYS:
            try:
                candidate = float(outputs[output_key])
            except (KeyError, TypeError, ValueError):
                continue
            if math.isfinite(candidate):
                mean_depth = candidate
                break
        if mean_depth is None:
            expected_keys = ", ".join(FIRE_MEAN_DEPTH_OUTPUT_KEYS)
            raise ValueError(
                f"{wildcards.sample}: {input.manifest} does not contain a "
                f"valid mean-depth output; expected one of: {expected_keys}"
            )

        if mean_depth < params.minimum_mean_depth:
            eligible = False
            reason = "below_minimum_mean_depth"
        else:
            eligible = True
            reason = "eligible"
        decision = {
            "sample": wildcards.sample,
            "mean_depth": mean_depth,
            "minimum_mean_depth": params.minimum_mean_depth,
            "eligible": eligible,
            "reason": reason,
        }
        with open(output.manifest, "w", encoding="utf-8") as handle:
            json.dump(decision, handle, indent=2)
            handle.write("\n")
        with open(output.report, "w", encoding="utf-8") as handle:
            handle.write(
                "sample\tmean_depth\tminimum_mean_depth\teligible\treason\n"
            )
            handle.write(
                f"{decision['sample']}\t{decision['mean_depth']:g}\t"
                f"{params.minimum_mean_depth:g}\t"
                f"{'yes' if decision['eligible'] else 'no'}\t"
                f"{decision['reason']}\n"
            )


rule fire_eligibility_report:
    input:
        rows=fire_eligibility_rows
    output:
        scratch=FIRE_ELIGIBILITY_REPORT,
        published=PUBLISHED_FIRE_ELIGIBILITY_REPORT
    run:
        header = (
            "sample\tmean_depth\tminimum_mean_depth\teligible\treason\n"
        )
        lines = [header]
        for row_path in input.rows:
            row_lines = Path(row_path).read_text(encoding="utf-8").splitlines()
            if not row_lines or f"{row_lines[0]}\n" != header:
                raise ValueError(f"invalid FIRE eligibility row: {row_path}")
            if len(row_lines) != 2:
                raise ValueError(
                    f"expected one FIRE eligibility decision in {row_path}"
                )
            lines.append(f"{row_lines[1]}\n")
        report = "".join(lines)
        for report_path in (output.scratch, output.published):
            destination = Path(report_path)
            destination.parent.mkdir(parents=True, exist_ok=True)
            temporary = destination.with_name(
                f".{destination.name}.{os.getpid()}.tmp"
            )
            temporary.write_text(report, encoding="utf-8")
            os.replace(temporary, destination)


rule copy_pbwdl:
    input:
        outputs=lambda wildcards: pbwgs_outputs_for_sample(wildcards.sample),
        tool=COPY_PBWDL,
        validation=reanalysis_validation_input
    output:
        bam=HAPLOTAGGED_BAM,
        manifest=PUBLISHED_WGS_MANIFEST
    threads: 8
    resources:
        mem_mb=32000,
        runtime=90,
        io_heavy=1
    log:
        "logs/copy_pbwdl/{sample}.log"
    benchmark:
        "benchmarks/copy_pbwdl/{sample}.tsv"
    params:
        output_directory="results/wgs/{sample}"
    shell:
        """
        mkdir -p {params.output_directory:q} "$(dirname {log:q})"
        {input.tool:q} \
            -i {input.outputs:q} \
            -o {params.output_directory:q} \
            -t {threads} -m > {log:q} 2>&1
        test -s {output.bam:q}
        test -s {output.manifest:q}
        """


rule fire_input_bam:
    input:
        bam=fire_source_bam
    output:
        bam=temp(FIRE_INPUT_BAM)
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10
    log:
        "logs/fire_input_bam/{sample}.log"
    benchmark:
        "benchmarks/fire_input_bam/{sample}.tsv"
    shell:
        """
        mkdir -p "$(dirname {output.bam:q})" "$(dirname {log:q})"
        source_bam="$(realpath {input.bam:q})"
        ln -sfn "$source_bam" {output.bam:q}
        test -s {output.bam:q}
        printf 'Using %s\n' "$source_bam" > {log:q}
        """


rule fire_input_bam_index:
    input:
        bam=FIRE_INPUT_BAM
    output:
        bai=temp(FIRE_INPUT_BAM_INDEX)
    threads: 8
    resources:
        mem_mb=16000,
        runtime=60,
        io_heavy=1
    log:
        "logs/fire_input_bam_index/{sample}.log"
    benchmark:
        "benchmarks/fire_input_bam_index/{sample}.tsv"
    conda:
        "envs/hifihla.yaml"
    shell:
        """
        mkdir -p "$(dirname {log:q})"
        samtools index -@ {threads} -o {output.bai:q} {input.bam:q} \
            > {log:q} 2>&1
        test -s {output.bai:q}
        """


rule haplotagged_bam_index:
    input:
        HAPLOTAGGED_BAM
    output:
        HAPLOTAGGED_BAM_INDEX
    threads: 8
    resources:
        mem_mb=16000,
        runtime=60,
        io_heavy=1
    log:
        "logs/haplotagged_bam_index/{sample}.log"
    benchmark:
        "benchmarks/haplotagged_bam_index/{sample}.tsv"
    conda:
        "envs/hifihla.yaml"
    shell:
        """
        mkdir -p "$(dirname {log:q})"
        samtools index -@ {threads} -o {output:q} {input:q} > {log:q} 2>&1
        test -s {output:q}
        """


include: "workflow/rules/fire.smk"


rule fire_bigwig:
    input:
        pileup=FIRE_PILEUP,
        index=FIRE_PILEUP_INDEX,
        chrom_sizes=CHROM_SIZES,
        shrink_tool=SHRINK_FIRE_BED
    output:
        FIRE_BIGWIG
    threads: 8
    resources:
        mem_mb=16000,
        runtime=60,
        io_heavy=1
    log:
        "logs/fire_bigwig/{sample}.log"
    benchmark:
        "benchmarks/fire_bigwig/{sample}.tsv"
    conda:
        "envs/intervals.yaml"
    shell:
        """
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        temporary_bedgraph="$(
            mktemp "${{TMPDIR:-/tmp}}/{wildcards.sample}.fire.XXXXXX"
        )"
        trap 'rm -f "$temporary_bedgraph"' EXIT
        (
            {input.shrink_tool:q} \
                -c {input.chrom_sizes:q} \
                {input.pileup:q} \
                "$temporary_bedgraph" \
                -t {threads}
            bedGraphToBigWig \
                "$temporary_bedgraph" \
                {input.chrom_sizes:q} \
                {output:q}
        ) > {log:q} 2>&1
        test -s {output:q}
        """


rule scatac_map:
    input:
        peaks=lambda wildcards: CCRE_PEAKS[wildcards.peak_set],
        pileup=FIRE_PILEUP,
        pileup_index=FIRE_PILEUP_INDEX,
        chrom_sizes=CHROM_SIZES
    output:
        temp("intermediate/fire/ccre/{peak_set}/{sample}.mapped.bed")
    threads: 1
    resources:
        mem_mb=16000,
        runtime=60
    params:
        sorted=int(SCATAC_SORTED)
    log:
        "logs/scatac_map/{peak_set}/{sample}.log"
    benchmark:
        "benchmarks/scatac_map/{peak_set}/{sample}.tsv"
    conda:
        "envs/intervals.yaml"
    shell:
        """
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        if [[ {params.sorted} -eq 1 ]]; then
            pileup_contigs="$(tabix --list-chroms {input.pileup:q})"
            bedtools map \
                -a {input.peaks:q} \
                -b <(
                    while read -r chrom _; do
                        if grep -Fqx "$chrom" <<< "$pileup_contigs"; then
                            tabix {input.pileup:q} "$chrom"
                        fi
                    done < {input.chrom_sizes:q}
                ) \
                -c 6,6,6 \
                -o max,median,sum \
                -sorted \
                -g {input.chrom_sizes:q} \
                > {output:q} 2> {log:q}
        else
            bedtools map \
                -a {input.peaks:q} \
                -b {input.pileup:q} \
                -c 6,6,6 \
                -o max,median,sum \
                > {output:q} 2> {log:q}
        fi
        test -s {output:q}
        """


rule scatac_split:
    input:
        mapped="intermediate/fire/ccre/{peak_set}/{sample}.mapped.bed"
    output:
        maximum=(
            "results/fire/ccre/{peak_set}/"
            "{sample}_scATAC_max_scores.bed"
        ),
        median=(
            "results/fire/ccre/{peak_set}/"
            "{sample}_scATAC_median_scores.bed"
        ),
        total=(
            "results/fire/ccre/{peak_set}/"
            "{sample}_scATAC_sum_scores.bed"
        )
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10
    benchmark:
        "benchmarks/scatac_split/{peak_set}/{sample}.tsv"
    conda:
        "envs/intervals.yaml"
    script:
        "workflow/scripts/split_scatac_map.py"


rule sawfish_discover:
    input:
        bam=HAPLOTAGGED_BAM,
        bai=HAPLOTAGGED_BAM + ".bai",
        reference=FASTA,
        reference_index=FASTA + ".fai",
        expected_cn=lambda wildcards: EXPECTED_CN[sex_of(wildcards.sample)],
        excluded=CNV_EXCLUDED
    output:
        directory("sawfish/{sample}/discover")
    threads: 24
    resources:
        mem_mb=64000,
        runtime=240
    log:
        "logs/sawfish_discover/{sample}.log"
    benchmark:
        "benchmarks/sawfish_discover/{sample}.tsv"
    conda:
        "envs/sawfish.yaml"
    shell:
        """
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        sawfish discover \
            --threads {threads} \
            --ref {input.reference:q} \
            --bam {input.bam:q} \
            --expected-cn {input.expected_cn:q} \
            --cnv-excluded-regions {input.excluded:q} \
            --output-dir {output:q} > {log:q} 2>&1
        """


rule sawfish_jointcall:
    input:
        discover="sawfish/{sample}/discover"
    output:
        directory("sawfish/{sample}/joint")
    threads: 24
    resources:
        mem_mb=64000,
        runtime=240
    log:
        "logs/sawfish_jointcall/{sample}.log"
    benchmark:
        "benchmarks/sawfish_jointcall/{sample}.tsv"
    conda:
        "envs/sawfish.yaml"
    shell:
        """
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        sawfish joint-call \
            --threads {threads} \
            --report-supporting-reads \
            --sample {input.discover:q} \
            --output-dir {output:q} > {log:q} 2>&1
        test -s {output:q}/genotyped.sv.vcf.gz
        test -s {output:q}/supporting_reads.json.gz
        """


rule publish_sawfish:
    input:
        joint="sawfish/{sample}/joint"
    output:
        vcf="results/sv/{sample}.sawfish.vcf.gz",
        index="results/sv/{sample}.sawfish.vcf.gz.tbi",
        reads="results/sv/{sample}.supporting_reads.json.gz"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30,
        io_heavy=1
    log:
        "logs/publish_sawfish/{sample}.log"
    benchmark:
        "benchmarks/publish_sawfish/{sample}.tsv"
    conda:
        "envs/sawfish.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.vcf:q})" "$(dirname {log:q})"
        (
            cp {input.joint:q}/genotyped.sv.vcf.gz {output.vcf:q}
            cp {input.joint:q}/supporting_reads.json.gz {output.reads:q}
            tabix -f -p vcf {output.vcf:q}
        ) > {log:q} 2>&1
        test -s {output.vcf:q}
        test -s {output.index:q}
        test -s {output.reads:q}
        """


rule svtopo:
    input:
        bam=HAPLOTAGGED_BAM,
        bai=HAPLOTAGGED_BAM + ".bai",
        vcf="results/sv/{sample}.sawfish.vcf.gz",
        vcf_index="results/sv/{sample}.sawfish.vcf.gz.tbi",
        reads="results/sv/{sample}.supporting_reads.json.gz",
        excluded=CNV_EXCLUDED
    output:
        directory("results/sv/{sample}.svtopo")
    threads: 1
    resources:
        mem_mb=6000,
        runtime=60
    params:
        prefix=lambda wildcards: svtopo_prefix(wildcards.sample)
    log:
        "logs/svtopo/{sample}.log"
    benchmark:
        "benchmarks/svtopo/{sample}.tsv"
    conda:
        "envs/svtopo.yaml"
    shell:
        """
        mkdir -p {output:q} "$(dirname {log:q})"
        svtopo \
            --bam {input.bam:q} \
            --svtopo-dir {output:q} \
            --prefix {params.prefix:q} \
            --variant-readnames {input.reads:q} \
            --vcf {input.vcf:q} \
            --exclude-regions {input.excluded:q} > {log:q} 2>&1
        find {output:q} -type f -size +0c -print -quit | grep -q .
        """


rule haplotagged_bam_pbindex:
    input:
        HAPLOTAGGED_BAM
    output:
        HAPLOTAGGED_BAM + ".pbi"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=60,
        io_heavy=1
    log:
        "logs/haplotagged_bam_pbindex/{sample}.log"
    benchmark:
        "benchmarks/haplotagged_bam_pbindex/{sample}.tsv"
    conda:
        "envs/pbtk.yaml"
    shell:
        """
        mkdir -p "$(dirname {log:q})"
        pbindex -j {threads} {input:q} > {log:q} 2>&1
        test -s {output:q}
        """


rule hla_fastq:
    input:
        bam=HAPLOTAGGED_BAM,
        pbi=HAPLOTAGGED_BAM + ".pbi"
    output:
        temp("intermediate/hla/{sample}/{sample}.fastq.gz")
    threads: 8
    resources:
        mem_mb=32000,
        runtime=180,
        io_heavy=1
    params:
        prefix="intermediate/hla/{sample}/{sample}"
    log:
        "logs/hla_fastq/{sample}.log"
    benchmark:
        "benchmarks/hla_fastq/{sample}.tsv"
    conda:
        "envs/pbtk.yaml"
    shell:
        """
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        bam2fastq \
            -j {threads} \
            {input.bam:q} \
            -o {params.prefix:q} > {log:q} 2>&1
        test -s {output:q}
        """


rule hla_hifiasm:
    input:
        "intermediate/hla/{sample}/{sample}.fastq.gz"
    output:
        hap1="intermediate/hla/{sample}/{sample}.asm.bp.hap1.p_ctg.gfa",
        hap2="intermediate/hla/{sample}/{sample}.asm.bp.hap2.p_ctg.gfa"
    threads: 32
    resources:
        mem_mb=384000,
        runtime=720
    params:
        prefix="intermediate/hla/{sample}/{sample}.asm"
    log:
        "logs/hla_hifiasm/{sample}.log"
    benchmark:
        "benchmarks/hla_hifiasm/{sample}.tsv"
    conda:
        "envs/hifiasm.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.hap1:q})" "$(dirname {log:q})"
        hifiasm \
            -o {params.prefix:q} \
            -t {threads} \
            {input:q} > {log:q} 2>&1
        test -s {output.hap1:q}
        test -s {output.hap2:q}
        """


rule hla_gfa_to_fasta:
    input:
        hap1="intermediate/hla/{sample}/{sample}.asm.bp.hap1.p_ctg.gfa",
        hap2="intermediate/hla/{sample}/{sample}.asm.bp.hap2.p_ctg.gfa"
    output:
        hap1="intermediate/hla/{sample}/{sample}.asm.bp.hap1.p_ctg.fa",
        hap2="intermediate/hla/{sample}/{sample}.asm.bp.hap2.p_ctg.fa"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=30
    log:
        "logs/hla_gfa_to_fasta/{sample}.log"
    benchmark:
        "benchmarks/hla_gfa_to_fasta/{sample}.tsv"
    shell:
        """
        mkdir -p "$(dirname {log:q})"
        (
            awk -F '\\t' '$1 == "S" {{print ">" $2 "\\n" $3}}' \
                {input.hap1:q} > {output.hap1:q}
            awk -F '\\t' '$1 == "S" {{print ">" $2 "\\n" $3}}' \
                {input.hap2:q} > {output.hap2:q}
        ) 2> {log:q}
        test -s {output.hap1:q}
        test -s {output.hap2:q}
        """


rule hla_align_contigs:
    input:
        hap1="intermediate/hla/{sample}/{sample}.asm.bp.hap1.p_ctg.fa",
        hap2="intermediate/hla/{sample}/{sample}.asm.bp.hap2.p_ctg.fa",
        reference=HLA_REF
    output:
        "intermediate/hla/{sample}/{sample}.GRCh38.asm_sorted.bam"
    threads: 20
    resources:
        mem_mb=64000,
        runtime=240
    params:
        align_threads=16,
        sort_threads=4,
        sort_memory="8G"
    log:
        "logs/hla_align_contigs/{sample}.log"
    benchmark:
        "benchmarks/hla_align_contigs/{sample}.tsv"
    conda:
        "envs/hifihla.yaml"
    shell:
        r"""
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        sort_prefix="${{TMPDIR:-/tmp}}/{wildcards.sample}.hla-sort"
        minimap2 \
            -t {params.align_threads} \
            -L --secondary=no --eqx -ax asm5 \
            -R '@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}' \
            {input.reference:q} \
            {input.hap1:q} {input.hap2:q} 2>> {log:q} |
        samtools sort \
            -@ {params.sort_threads} \
            -m {params.sort_memory} \
            -T "$sort_prefix" \
            -o {output:q} 2>> {log:q}
        test -s {output:q}
        """


rule hla_contig_bam_index:
    input:
        "intermediate/hla/{sample}/{sample}.GRCh38.asm_sorted.bam"
    output:
        "intermediate/hla/{sample}/{sample}.GRCh38.asm_sorted.bam.bai"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=60,
        io_heavy=1
    log:
        "logs/hla_contig_bam_index/{sample}.log"
    benchmark:
        "benchmarks/hla_contig_bam_index/{sample}.tsv"
    conda:
        "envs/hifihla.yaml"
    shell:
        """
        mkdir -p "$(dirname {log:q})"
        samtools index -@ {threads} -o {output:q} {input:q} > {log:q} 2>&1
        test -s {output:q}
        """


if HLA_MODE == "reads":

    rule hla_call:
        input:
            bam=HAPLOTAGGED_BAM,
            bai=HAPLOTAGGED_BAM + ".bai"
        output:
            json=HLA_JSON,
            report=HLA_REPORT,
            summary=HLA_SUMMARY
        threads: 8
        resources:
            mem_mb=32000,
            runtime=120
        params:
            prefix="results/hla/{sample}"
        log:
            "logs/hla_call/{sample}.log"
        benchmark:
            "benchmarks/hla_call/{sample}.tsv"
        conda:
            "envs/hifihla.yaml"
        shell:
            """
            mkdir -p "$(dirname {output.json:q})" "$(dirname {log:q})"
            hifihla call-reads \
                --preset wgs \
                -j {threads} \
                -a {input.bam:q} \
                -o {params.prefix:q} > {log:q} 2>&1
            test -s {output.json:q}
            test -s {output.report:q}
            test -s {output.summary:q}
            """

else:

    rule hla_call:
        input:
            bam="intermediate/hla/{sample}/{sample}.GRCh38.asm_sorted.bam",
            bai="intermediate/hla/{sample}/{sample}.GRCh38.asm_sorted.bam.bai",
            hap1="intermediate/hla/{sample}/{sample}.asm.bp.hap1.p_ctg.fa",
            hap2="intermediate/hla/{sample}/{sample}.asm.bp.hap2.p_ctg.fa"
        output:
            json=HLA_JSON,
            report=HLA_REPORT,
            summary=HLA_SUMMARY
        threads: 8
        resources:
            mem_mb=32000,
            runtime=120
        params:
            prefix="results/hla/{sample}"
        log:
            "logs/hla_call/{sample}.log"
        benchmark:
            "benchmarks/hla_call/{sample}.tsv"
        conda:
            "envs/hifihla.yaml"
        shell:
            """
            mkdir -p "$(dirname {output.json:q})" "$(dirname {log:q})"
            hifihla call-contigs \
                --abam {input.bam:q} \
                --hap1 {input.hap1:q} \
                --hap2 {input.hap2:q} \
                --threads {threads} \
                --out_prefix {params.prefix:q} > {log:q} 2>&1
            test -s {output.json:q}
            test -s {output.report:q}
            test -s {output.summary:q}
            """
