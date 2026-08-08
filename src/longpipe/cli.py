"""Manage durable Longpipe runs backed by temporary scratch storage."""

from __future__ import annotations

import argparse
import csv
import fcntl
import json
import os
import re
import shutil
import signal
import socket
import subprocess
import sys
import sysconfig
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterator, Sequence

import yaml

from . import __version__


DEFAULT_PROJECT_ROOT = "/fh/fast/furlan_s/project"
DEFAULT_SCRATCH_ROOT = "/hpc/temp/furlan_s"
RUN_NAME_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


class LongpipeError(RuntimeError):
    """An expected command-line error."""


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _path(value: str | Path, base: Path | None = None) -> Path:
    expanded = Path(os.path.expandvars(os.path.expanduser(str(value))))
    if not expanded.is_absolute() and base is not None:
        expanded = base / expanded
    return expanded.resolve()


def _workflow_root() -> Path:
    override = os.environ.get("LONGPIPE_WORKFLOW_ROOT")
    candidates = []
    if override:
        candidates.append(_path(override))

    # Source checkout, followed by the wheel/data-files installation location.
    candidates.extend(
        [
            Path(__file__).resolve().parents[2],
            Path(sysconfig.get_path("data")) / "share" / "longpipe",
        ]
    )
    for candidate in candidates:
        if (
            (candidate / "Snakefile").is_file()
            and (candidate / "config" / "config.yaml").is_file()
            and (candidate / "profiles" / "slurm" / "config.yaml").is_file()
        ):
            return candidate
    checked = ", ".join(str(candidate) for candidate in candidates)
    raise LongpipeError(
        f"could not locate Longpipe workflow resources; checked {checked}"
    )


def _validate_run_name(run_name: str) -> None:
    if not RUN_NAME_PATTERN.fullmatch(run_name):
        raise LongpipeError(
            "run name must start with a letter or digit and contain only "
            "letters, digits, '.', '_', and '-'"
        )


def _validate_working_directory(project_root_value: str) -> tuple[Path, Path]:
    project_root = _path(project_root_value)
    working_directory = Path.cwd().resolve()
    if not project_root.is_dir():
        raise LongpipeError(f"configured project root does not exist: {project_root}")
    if (
        working_directory != project_root
        and project_root not in working_directory.parents
    ):
        raise LongpipeError(
            f"working directory must be within project root {project_root}; "
            f"got {working_directory}"
        )
    return project_root, working_directory


def _run_directory(project_root_value: str, run_name: str) -> Path:
    _validate_run_name(run_name)
    project_root, working_directory = _validate_working_directory(project_root_value)
    run_directory = working_directory / run_name
    if run_directory.exists() or run_directory.is_symlink():
        resolved = run_directory.resolve()
        if resolved != project_root and project_root not in resolved.parents:
            raise LongpipeError(
                f"run directory resolves outside project root {project_root}: "
                f"{run_directory}"
            )
        return resolved
    return run_directory


def _scratch_directory(scratch_root_value: str, run_name: str) -> tuple[Path, Path]:
    scratch_root = _path(scratch_root_value)
    scratch_directory = scratch_root / run_name
    if scratch_directory.parent != scratch_root:
        raise LongpipeError(f"invalid scratch path for run {run_name!r}")
    return scratch_root, scratch_directory


def _atomic_write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        with temporary.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


@contextmanager
def _run_lock(scratch_root: Path, run_name: str, run_directory: Path) -> Iterator[None]:
    lock_directory = scratch_root / ".longpipe-locks"
    lock_directory.mkdir(parents=True, exist_ok=True)
    lock_path = lock_directory / f"{run_name}.lock"
    with lock_path.open("a+", encoding="utf-8") as lock_handle:
        try:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise LongpipeError(
                f"run {run_name!r} is already active; lock: {lock_path}"
            ) from error

        lock_handle.seek(0)
        lock_handle.truncate()
        json.dump(
            {
                "hostname": socket.gethostname(),
                "pid": os.getpid(),
                "run_directory": str(run_directory),
                "run_name": run_name,
                "started_at": _now(),
            },
            lock_handle,
            indent=2,
            sort_keys=True,
        )
        lock_handle.write("\n")
        lock_handle.flush()
        os.fsync(lock_handle.fileno())
        try:
            yield
        finally:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)


def _owner_path(scratch_directory: Path) -> Path:
    return scratch_directory / ".longpipe" / "owner.json"


def _read_owner(scratch_directory: Path) -> dict[str, object]:
    owner_path = _owner_path(scratch_directory)
    if not owner_path.is_file():
        raise LongpipeError(
            f"scratch directory exists without Longpipe ownership metadata: "
            f"{scratch_directory}"
        )
    try:
        with owner_path.open(encoding="utf-8") as handle:
            owner = json.load(handle)
    except (OSError, json.JSONDecodeError) as error:
        raise LongpipeError(
            f"invalid scratch ownership metadata: {owner_path}"
        ) from error
    if not isinstance(owner, dict):
        raise LongpipeError(f"invalid scratch ownership metadata: {owner_path}")
    return owner


def _verify_owner(
    scratch_directory: Path, run_name: str, run_directory: Path
) -> dict[str, object]:
    if scratch_directory.is_symlink():
        raise LongpipeError(
            f"refusing to use symlinked scratch path: {scratch_directory}"
        )
    owner = _read_owner(scratch_directory)
    expected_directory = str(run_directory)
    if (
        owner.get("run_name") != run_name
        or owner.get("run_directory") != expected_directory
    ):
        raise LongpipeError(
            f"scratch directory {scratch_directory} belongs to "
            f"{owner.get('run_directory', 'an unknown run')}, not {run_directory}"
        )
    return owner


def _create_scratch(
    scratch_directory: Path, run_name: str, run_directory: Path
) -> None:
    if scratch_directory.exists():
        _verify_owner(scratch_directory, run_name, run_directory)
        return
    scratch_directory.mkdir(parents=True)
    _atomic_write_json(
        _owner_path(scratch_directory),
        {
            "created_at": _now(),
            "run_directory": str(run_directory),
            "run_name": run_name,
        },
    )


def _remove_scratch(scratch_root: Path, scratch_directory: Path) -> None:
    if scratch_directory.is_symlink():
        raise LongpipeError(
            f"refusing to remove symlinked scratch path: {scratch_directory}"
        )
    if scratch_directory.parent != scratch_root:
        raise LongpipeError(
            f"refusing to remove unsafe scratch path: {scratch_directory}"
        )
    shutil.rmtree(scratch_directory)


def _load_run_config(config_path: Path, run_directory: Path) -> tuple[dict, Path]:
    try:
        with config_path.open(encoding="utf-8") as handle:
            config = yaml.safe_load(handle)
    except OSError as error:
        raise LongpipeError(f"cannot read configuration: {config_path}") from error
    except yaml.YAMLError as error:
        raise LongpipeError(f"invalid YAML in {config_path}: {error}") from error

    if not isinstance(config, dict):
        raise LongpipeError(f"configuration must contain a YAML mapping: {config_path}")
    reserved = sorted({"name", "outdir"} & set(config))
    if reserved:
        raise LongpipeError(
            f"{config_path}: {', '.join(reserved)} are managed by the CLI and "
            "must not be set"
        )

    runsheet_value = config.get("runsheet")
    if not isinstance(runsheet_value, str) or not runsheet_value.strip():
        raise LongpipeError(f"{config_path}: runsheet must be a non-empty path")
    runsheet_path = _path(runsheet_value, run_directory)
    if not runsheet_path.is_file():
        raise LongpipeError(f"runsheet does not exist: {runsheet_path}")
    return config, runsheet_path


def _invocation_id() -> str:
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    return f"{timestamp}-{os.getpid()}-{uuid.uuid4().hex[:8]}"


def _snapshot_inputs(
    run_directory: Path,
    scratch_directory: Path,
    invocation_id: str,
) -> tuple[Path, Path, Path]:
    config_path = run_directory / "config.yaml"
    _, runsheet_path = _load_run_config(config_path, run_directory)
    invocation_directory = (
        scratch_directory / ".longpipe" / "invocations" / invocation_id
    )
    invocation_directory.mkdir(parents=True)
    config_snapshot = invocation_directory / "config.yaml"
    runsheet_snapshot = invocation_directory / "runsheet.csv"
    shutil.copy2(config_path, config_snapshot)
    shutil.copy2(runsheet_path, runsheet_snapshot)

    durable_invocation = (
        run_directory / "run-metadata" / "invocations" / invocation_id
    )
    durable_invocation.mkdir(parents=True)
    shutil.copy2(config_snapshot, durable_invocation / "config.yaml")
    shutil.copy2(runsheet_snapshot, durable_invocation / "runsheet.csv")
    return config_snapshot, runsheet_snapshot, durable_invocation


def _record_status(
    status: dict[str, object],
    scratch_directory: Path,
    run_directory: Path,
    durable_invocation: Path,
) -> None:
    _atomic_write_json(scratch_directory / ".longpipe" / "status.json", status)
    _atomic_write_json(durable_invocation / "status.json", status)
    _atomic_write_json(run_directory / "run-metadata" / "status.json", status)


def _replace_tree(
    source: Path,
    destination: Path,
    required: bool = False,
    prepare: Callable[[Path], None] | None = None,
) -> None:
    if required and not source.is_dir():
        raise LongpipeError(f"expected workflow output directory is missing: {source}")

    def build(stage: Path) -> None:
        if source.is_dir():
            shutil.copytree(source, stage)
        else:
            stage.mkdir()
        if prepare is not None:
            prepare(stage)

    _replace_built_tree(destination, build)


def _replace_built_tree(destination: Path, build: Callable[[Path], None]) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    stage = destination.parent / f".{destination.name}.stage-{uuid.uuid4().hex}"
    backup = destination.parent / f".{destination.name}.backup-{uuid.uuid4().hex}"
    try:
        build(stage)

        had_destination = destination.exists() or destination.is_symlink()
        if had_destination:
            os.replace(destination, backup)
        try:
            os.replace(stage, destination)
        except BaseException:
            if had_destination and backup.exists():
                os.replace(backup, destination)
            raise
        if backup.exists():
            if backup.is_dir() and not backup.is_symlink():
                shutil.rmtree(backup)
            else:
                backup.unlink()
    finally:
        if stage.exists():
            shutil.rmtree(stage)
        if backup.exists():
            if destination.exists():
                if backup.is_dir() and not backup.is_symlink():
                    shutil.rmtree(backup)
                else:
                    backup.unlink()
            else:
                os.replace(backup, destination)


def _atomic_copy(source: Path, destination: Path, required: bool = False) -> None:
    if required and not source.is_file():
        raise LongpipeError(f"expected workflow output file is missing: {source}")
    if not source.is_file():
        return

    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(
        f".{destination.name}.{uuid.uuid4().hex}.tmp"
    )
    try:
        shutil.copy2(source, temporary)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def _sample_labels(runsheet: Path) -> list[str]:
    try:
        with runsheet.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
    except OSError as error:
        raise LongpipeError(f"could not read runsheet {runsheet}: {error}") from error

    labels: list[str] = []
    for row in rows:
        label = (
            row.get("lab_name", "").strip()
            or row.get("core_name", "").strip()
        )
        if not RUN_NAME_PATTERN.fullmatch(label):
            raise LongpipeError(f"invalid sample label in {runsheet}: {label!r}")
        labels.append(label)
    if not labels:
        raise LongpipeError(f"runsheet contains no samples: {runsheet}")
    duplicates = sorted({label for label in labels if labels.count(label) > 1})
    if duplicates:
        raise LongpipeError(
            f"duplicate sample labels in {runsheet}: {', '.join(duplicates)}"
        )
    return labels


def _read_fire_eligibility(report: Path, samples: set[str]) -> set[str]:
    if not report.is_file():
        raise LongpipeError(f"expected workflow output file is missing: {report}")
    try:
        with report.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except OSError as error:
        raise LongpipeError(
            f"could not read FIRE eligibility report {report}: {error}"
        ) from error

    eligible: set[str] = set()
    for row in rows:
        sample = row.get("sample", "").strip()
        decision = row.get("eligible", "").strip().lower()
        if sample not in samples:
            raise LongpipeError(
                f"FIRE eligibility report contains unknown sample {sample!r}: {report}"
            )
        if decision not in {"yes", "no"}:
            raise LongpipeError(
                f"FIRE eligibility report has invalid decision for {sample!r}: "
                f"{decision!r}"
            )
        if decision == "yes":
            eligible.add(sample)
    return eligible


def _sample_from_filename(
    filename: str, known_samples: set[str], delimiter: str
) -> str | None:
    matches = [
        sample
        for sample in known_samples
        if filename.startswith(f"{sample}{delimiter}")
    ]
    return max(matches, key=len) if matches else None


def _copy_output(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if source.is_dir():
        shutil.copytree(source, destination)
    else:
        shutil.copy2(source, destination)


def _copy_required_output(source: Path, destination: Path) -> None:
    if not source.exists():
        raise LongpipeError(f"expected workflow output is missing: {source}")
    _copy_output(source, destination)


def _rewrite_published_wgs_manifest(
    manifest: Path, scratch_sample: Path, durable_sample: Path
) -> None:
    scratch_root = str(scratch_sample)
    durable_root = str(durable_sample)
    scratch_prefix = scratch_root + os.sep

    def rewrite(value: object) -> object:
        if isinstance(value, dict):
            return {key: rewrite(item) for key, item in value.items()}
        if isinstance(value, list):
            return [rewrite(item) for item in value]
        if isinstance(value, str):
            if value == scratch_root:
                return durable_root
            if value.startswith(scratch_prefix):
                return durable_root + os.sep + value[len(scratch_prefix) :]
        return value

    try:
        payload = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise LongpipeError(
            f"could not read WGS output manifest {manifest}: {error}"
        ) from error
    if not isinstance(payload, dict):
        raise LongpipeError(f"WGS output manifest is not a JSON object: {manifest}")
    rewritten = rewrite(payload)
    assert isinstance(rewritten, dict)
    _atomic_write_json(manifest, rewritten)


def _copy_fire_sample_outputs(
    scratch_fire: Path,
    destination: Path,
    sample: str,
    peak_sets: set[str] | None,
    fire_version: str,
) -> None:
    pileup_name = f"{sample}-fire-{fire_version}-pileup.bed.gz"
    _copy_required_output(
        scratch_fire / "pileups" / pileup_name,
        destination / "pileups" / pileup_name,
    )
    _copy_required_output(
        scratch_fire / "pileups" / f"{pileup_name}.tbi",
        destination / "pileups" / f"{pileup_name}.tbi",
    )
    bigwig_name = f"{sample}-fire-{fire_version}-pileup_shrunk.bw"
    _copy_required_output(
        scratch_fire / "bigwigs" / bigwig_name,
        destination / "bigwigs" / bigwig_name,
    )
    track_names = (
        f"{sample}-fire-{fire_version}-H1-percent-accessible.bw",
        f"{sample}-fire-{fire_version}-H2-percent-accessible.bw",
        (
            f"{sample}-fire-{fire_version}-"
            "H1-minus-H2-percent-accessible.bw"
        ),
    )
    for track_name in track_names:
        _copy_required_output(
            scratch_fire / "haplotypes" / "tracks" / track_name,
            destination / "haplotypes" / "tracks" / track_name,
        )
    differences_name = (
        f"{sample}-fire-{fire_version}-haplotype-differences.bed.gz"
    )
    difference_files = (
        differences_name,
        f"{differences_name}.tbi",
        (
            f"{sample}-fire-{fire_version}-"
            "haplotype-differences-summary.tsv"
        ),
    )
    for difference_name in difference_files:
        _copy_required_output(
            scratch_fire / "haplotypes" / "differences" / difference_name,
            destination / "haplotypes" / "differences" / difference_name,
        )
    if peak_sets is None:
        return
    for peak_set in peak_sets:
        for statistic in ("max", "median", "sum"):
            ccre_name = f"{sample}_scATAC_{statistic}_scores.bed"
            _copy_required_output(
                scratch_fire / "ccre" / peak_set / ccre_name,
                destination / "ccre" / peak_set / ccre_name,
            )


def _build_sample_results(
    stage: Path,
    scratch_results: Path,
    durable_samples: Path,
    sample_names: list[str],
    peak_sets: set[str],
    fire_version: str,
) -> None:
    expected_result_directories = {"fire", "hla", "sv", "wgs"}
    unexpected_result_entries = sorted(
        entry
        for entry in scratch_results.iterdir()
        if entry.name not in expected_result_directories
    )
    if unexpected_result_entries:
        raise LongpipeError(
            f"unexpected workflow result: {unexpected_result_entries[0]}"
        )

    scratch_wgs = scratch_results / "wgs"
    if not scratch_wgs.is_dir():
        raise LongpipeError(
            f"expected workflow output directory is missing: {scratch_wgs}"
        )

    scratch_sample_directories = {
        entry.name for entry in scratch_wgs.iterdir() if entry.is_dir()
    }
    unexpected_wgs_entries = [
        entry for entry in scratch_wgs.iterdir() if not entry.is_dir()
    ]
    if unexpected_wgs_entries:
        raise LongpipeError(
            f"unexpected entry in WGS result directory: {unexpected_wgs_entries[0]}"
        )

    current_samples = set(sample_names)
    missing_samples = sorted(current_samples - scratch_sample_directories)
    if missing_samples:
        raise LongpipeError(
            "missing WGS results for samples: " + ", ".join(missing_samples)
        )

    eligibility_report = scratch_results / "fire" / "eligibility.tsv"
    eligible_fire_samples = _read_fire_eligibility(
        eligibility_report, current_samples
    )
    stage.mkdir()

    for sample in sample_names:
        source = scratch_wgs / sample
        destination = stage / sample / "pbWGS"
        shutil.copytree(source, destination)
        manifest = destination / "outputs.json"
        if not manifest.is_file():
            raise LongpipeError(f"WGS output manifest is missing: {source}")
        _rewrite_published_wgs_manifest(
            manifest,
            source,
            durable_samples / sample / "pbWGS",
        )

    scratch_fire = scratch_results / "fire"
    expected_fire_entries = {
        "bigwigs",
        "ccre",
        "eligibility.tsv",
        "haplotypes",
        "pileups",
    }
    unexpected_fire_entries = sorted(
        entry
        for entry in scratch_fire.iterdir()
        if entry.name not in expected_fire_entries
    )
    if unexpected_fire_entries:
        raise LongpipeError(
            f"unexpected FIRE result: {unexpected_fire_entries[0]}"
        )
    for sample in eligible_fire_samples:
        _copy_fire_sample_outputs(
            scratch_fire,
            stage / sample / "fire",
            sample,
            peak_sets,
            fire_version,
        )

    scratch_sv = scratch_results / "sv"
    if not scratch_sv.is_dir():
        raise LongpipeError(
            f"expected workflow output directory is missing: {scratch_sv}"
        )
    for source in sorted(scratch_sv.iterdir()):
        sample = _sample_from_filename(
            source.name, scratch_sample_directories, "."
        )
        if sample is None:
            raise LongpipeError(f"could not assign SV output to a sample: {source}")
        if sample not in current_samples:
            continue
        if source.name == f"{sample}.svtopo" and source.is_dir():
            destination = stage / sample / "svtopo"
        else:
            destination = stage / sample / "sawfish" / source.name
        _copy_output(source, destination)

    scratch_hla = scratch_results / "hla"
    if not scratch_hla.is_dir():
        raise LongpipeError(
            f"expected workflow output directory is missing: {scratch_hla}"
        )
    for source in sorted(scratch_hla.iterdir()):
        sample = _sample_from_filename(
            source.name, scratch_sample_directories, "_"
        )
        if sample is None:
            raise LongpipeError(f"could not assign HLA output to a sample: {source}")
        if sample in current_samples:
            _copy_output(source, stage / sample / "hla" / source.name)

    for sample in sample_names:
        expected = (
            stage / sample / "sawfish" / f"{sample}.sawfish.vcf.gz",
            stage / sample / "sawfish" / f"{sample}.sawfish.vcf.gz.tbi",
            stage / sample / "sawfish" / f"{sample}.supporting_reads.json.gz",
            stage / sample / "svtopo",
            stage / sample / "hla" / f"{sample}_hifihla_report.json",
            stage / sample / "hla" / f"{sample}_hifihla_report.tsv",
            stage / sample / "hla" / f"{sample}_hifihla_summary.tsv",
        )
        missing = [str(path) for path in expected if not path.exists()]
        if missing:
            raise LongpipeError(
                f"missing published outputs for {sample}: {', '.join(missing)}"
            )


def _publish_sample_results(
    scratch_results: Path,
    run_directory: Path,
    runsheet_snapshot: Path,
    config: dict,
) -> None:
    if not scratch_results.is_dir():
        raise LongpipeError(
            f"expected workflow output directory is missing: {scratch_results}"
        )
    sample_names = _sample_labels(runsheet_snapshot)
    ccre_peaks = config.get("ccre_peaks", {})
    if not isinstance(ccre_peaks, dict):
        raise LongpipeError("ccre_peaks must be a mapping")
    fire_version = config.get("fire_version")
    if not isinstance(fire_version, str) or not fire_version:
        raise LongpipeError("fire_version must be a non-empty string")
    durable_samples = run_directory / "samples"

    _replace_built_tree(
        durable_samples,
        lambda stage: _build_sample_results(
            stage,
            scratch_results,
            durable_samples,
            sample_names,
            set(ccre_peaks),
            fire_version,
        ),
    )
    _atomic_copy(
        scratch_results / "fire" / "eligibility.tsv",
        run_directory / "run-metadata" / "fire-eligibility.tsv",
        required=True,
    )

    legacy_results = run_directory / "results"
    if legacy_results.is_symlink() or legacy_results.is_file():
        legacy_results.unlink()
    elif legacy_results.is_dir():
        shutil.rmtree(legacy_results)


def _publish_fire_only_results(
    scratch_results: Path,
    run_directory: Path,
    runsheet_snapshot: Path,
    config: dict,
) -> None:
    if not scratch_results.is_dir():
        raise LongpipeError(
            f"expected workflow output directory is missing: {scratch_results}"
        )
    sample_names = _sample_labels(runsheet_snapshot)
    fire_version = config.get("fire_version")
    if not isinstance(fire_version, str) or not fire_version:
        raise LongpipeError("fire_version must be a non-empty string")
    scratch_fire = scratch_results / "fire"
    if not scratch_fire.is_dir():
        raise LongpipeError(
            f"expected workflow output directory is missing: {scratch_fire}"
        )
    allowed_fire_entries = {
        "bigwigs",
        "ccre",
        "eligibility.tsv",
        "haplotypes",
        "pileups",
    }
    unexpected_fire_entries = sorted(
        entry
        for entry in scratch_fire.iterdir()
        if entry.name not in allowed_fire_entries
    )
    if unexpected_fire_entries:
        raise LongpipeError(
            f"unexpected FIRE result: {unexpected_fire_entries[0]}"
        )
    eligible_samples = _read_fire_eligibility(
        scratch_fire / "eligibility.tsv",
        set(sample_names),
    )
    durable_samples = run_directory / "samples"
    for sample in sample_names:
        destination = durable_samples / sample / "fire"
        if sample not in eligible_samples:
            if destination.is_symlink() or destination.is_file():
                destination.unlink()
            elif destination.is_dir():
                shutil.rmtree(destination)
            continue

        existing_ccre = destination / "ccre"

        def build(
            stage: Path,
            sample_name: str = sample,
            preserved_ccre: Path = existing_ccre,
        ) -> None:
            stage.mkdir()
            if preserved_ccre.is_dir():
                shutil.copytree(preserved_ccre, stage / "ccre")
            _copy_fire_sample_outputs(
                scratch_fire,
                stage,
                sample_name,
                None,
                fire_version,
            )

        _replace_built_tree(destination, build)

    _atomic_copy(
        scratch_fire / "eligibility.tsv",
        run_directory / "run-metadata" / "fire-eligibility.tsv",
        required=True,
    )


def _publish_snakemake_metadata(
    scratch_directory: Path, durable_invocation: Path
) -> None:
    source_root = scratch_directory / ".snakemake"
    destination_root = durable_invocation / "snakemake"
    destination_root.mkdir(exist_ok=True)
    for name in ("log", "metadata", "incomplete"):
        source = source_root / name
        if source.is_dir():
            shutil.copytree(source, destination_root / name, dirs_exist_ok=True)


def _publish_diagnostics(
    scratch_directory: Path, run_directory: Path, durable_invocation: Path
) -> None:
    _replace_tree(scratch_directory / "logs", run_directory / "logs")
    _replace_tree(scratch_directory / "benchmarks", run_directory / "benchmarks")
    _publish_snakemake_metadata(scratch_directory, durable_invocation)


def _snakemake_command(
    executable: str,
    workflow_root: Path,
    config_snapshot: Path,
    targets: Sequence[str] = (),
) -> list[str]:
    return [
        executable,
        *targets,
        "--snakefile",
        str(workflow_root / "Snakefile"),
        "--profile",
        str(workflow_root / "profiles" / "slurm"),
        "--configfile",
        str(config_snapshot),
    ]


def _resolve_snakemake_executable(explicit: str | None) -> str:
    if explicit:
        return explicit

    environment_override = os.environ.get("LONGPIPE_SNAKEMAKE")
    if environment_override:
        return environment_override

    sibling = Path(sys.executable).with_name("snakemake")
    if sibling.is_file() and os.access(sibling, os.X_OK):
        return str(sibling)

    executable = shutil.which("snakemake")
    if executable:
        return executable

    raise LongpipeError(
        "could not locate Snakemake; reinstall Longpipe with Pipx or use "
        "--snakemake-executable"
    )


def _stream_command(
    command: Sequence[str],
    environment: dict[str, str],
    working_directory: Path,
    log_path: Path,
) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as log_handle:
        log_handle.write(f"\n[{_now()}] command: {' '.join(command)}\n")
        log_handle.flush()
        try:
            process = subprocess.Popen(
                command,
                cwd=working_directory,
                env=environment,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                errors="replace",
                bufsize=1,
            )
        except OSError as error:
            raise LongpipeError(f"could not start {command[0]}: {error}") from error

        assert process.stdout is not None
        try:
            for line in process.stdout:
                sys.stdout.write(line)
                sys.stdout.flush()
                log_handle.write(line)
                log_handle.flush()
            return process.wait()
        except KeyboardInterrupt:
            process.send_signal(signal.SIGINT)
            try:
                process.wait(timeout=30)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
            return 130


def _setup(args: argparse.Namespace) -> int:
    run_directory = _run_directory(args.project_root, args.run_name)
    if run_directory.exists():
        raise LongpipeError(f"run directory already exists: {run_directory}")

    workflow_root = _workflow_root()
    config_template = workflow_root / "config" / "config.yaml"
    runsheet_template = workflow_root / "config" / "runsheet.csv"
    try:
        run_directory.mkdir()
        shutil.copy2(config_template, run_directory / "config.yaml")
        shutil.copy2(runsheet_template, run_directory / "runsheet.csv")
    except BaseException:
        if run_directory.exists():
            shutil.rmtree(run_directory)
        raise

    print(f"Created run {args.run_name!r} in {run_directory}")
    print(f"Edit {run_directory / 'config.yaml'} and {run_directory / 'runsheet.csv'}")
    return 0


def _run(args: argparse.Namespace) -> int:
    run_directory = _run_directory(args.project_root, args.run_name)
    if not run_directory.is_dir():
        raise LongpipeError(f"run directory does not exist: {run_directory}")
    _load_run_config(run_directory / "config.yaml", run_directory)

    scratch_root, scratch_directory = _scratch_directory(
        args.scratch_root, args.run_name
    )
    workflow_root = _workflow_root()
    snakemake_executable = _resolve_snakemake_executable(
        args.snakemake_executable
    )

    with _run_lock(scratch_root, args.run_name, run_directory):
        if scratch_directory.exists():
            _verify_owner(scratch_directory, args.run_name, run_directory)
            if args.restart:
                print(f"Removing existing scratch run {scratch_directory}")
                _remove_scratch(scratch_root, scratch_directory)
        _create_scratch(scratch_directory, args.run_name, run_directory)

        invocation_id = _invocation_id()
        config_snapshot, runsheet_snapshot, durable_invocation = _snapshot_inputs(
            run_directory, scratch_directory, invocation_id
        )
        config, _ = _load_run_config(config_snapshot, config_snapshot.parent)
        run_mode = "fire_only" if args.only_fire else "full"
        targets = ("fire_only",) if args.only_fire else ()
        command = _snakemake_command(
            snakemake_executable,
            workflow_root,
            config_snapshot,
            targets,
        )
        status: dict[str, object] = {
            "command": command,
            "hostname": socket.gethostname(),
            "invocation_id": invocation_id,
            "pid": os.getpid(),
            "run_directory": str(run_directory),
            "mode": run_mode,
            "run_name": args.run_name,
            "scratch_directory": str(scratch_directory),
            "started_at": _now(),
            "status": "running",
        }
        _record_status(status, scratch_directory, run_directory, durable_invocation)

        environment = os.environ.copy()
        environment.update(
            {
                "LONGPIPE_INVOCATION_ID": invocation_id,
                "LONGPIPE_RUN_DIR": str(run_directory),
                "LONGPIPE_RUN_NAME": args.run_name,
                "LONGPIPE_RUNSHEET": str(runsheet_snapshot),
                "LONGPIPE_SCRATCH_DIR": str(scratch_directory),
            }
        )

        mode_label = "FIRE only" if args.only_fire else "full workflow"
        print(f"Running {mode_label} for {args.run_name!r} in {scratch_directory}")
        try:
            return_code = _stream_command(
                command,
                environment,
                run_directory,
                scratch_directory / "logs" / "longpipe" / "driver.log",
            )
        except LongpipeError as error:
            return_code = 1
            status["error"] = str(error)

        status["finished_at"] = _now()
        status["return_code"] = return_code
        status["status"] = "publishing" if return_code == 0 else "failed"
        _record_status(status, scratch_directory, run_directory, durable_invocation)

        try:
            _publish_diagnostics(
                scratch_directory, run_directory, durable_invocation
            )
            if return_code == 0:
                scratch_results = scratch_directory / "results"
                if args.only_fire:
                    _publish_fire_only_results(
                        scratch_results,
                        run_directory,
                        runsheet_snapshot,
                        config,
                    )
                else:
                    _publish_sample_results(
                        scratch_results,
                        run_directory,
                        runsheet_snapshot,
                        config,
                    )
        except BaseException as error:
            status["error"] = f"publication failed: {error}"
            status["status"] = "publication_failed"
            status["finished_at"] = _now()
            _record_status(status, scratch_directory, run_directory, durable_invocation)
            raise LongpipeError(str(status["error"])) from error

        status["status"] = "success" if return_code == 0 else "failed"
        status["finished_at"] = _now()
        _record_status(status, scratch_directory, run_directory, durable_invocation)
        if return_code == 0:
            completion = "FIRE-only run" if args.only_fire else "Run"
            print(
                f"{completion} completed; results published to "
                f"{run_directory / 'samples'}"
            )
        else:
            print(
                f"Run failed with exit code {return_code}; diagnostics published to "
                f"{run_directory}",
                file=sys.stderr,
            )
        return return_code


def _clean(args: argparse.Namespace) -> int:
    run_directory = _run_directory(args.project_root, args.run_name)
    if not run_directory.is_dir():
        raise LongpipeError(f"run directory does not exist: {run_directory}")
    scratch_root, scratch_directory = _scratch_directory(
        args.scratch_root, args.run_name
    )
    if not scratch_directory.exists():
        print(f"No scratch data exists for {args.run_name!r}")
        return 0

    with _run_lock(scratch_root, args.run_name, run_directory):
        _verify_owner(scratch_directory, args.run_name, run_directory)
        _remove_scratch(scratch_root, scratch_directory)
    print(f"Removed scratch data for {args.run_name!r}: {scratch_directory}")
    return 0


def _common_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "--project-root",
        default=os.environ.get("LONGPIPE_PROJECT_ROOT", DEFAULT_PROJECT_ROOT),
        help=(
            "allowed durable project root "
            "(default: LONGPIPE_PROJECT_ROOT or %(default)s)"
        ),
    )
    return parser


def _scratch_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "--scratch-root",
        default=os.environ.get("LONGPIPE_SCRATCH_ROOT", DEFAULT_SCRATCH_ROOT),
        help=(
            "temporary workflow root "
            "(default: LONGPIPE_SCRATCH_ROOT or %(default)s)"
        ),
    )
    return parser


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="longpipe",
        description="Run Longpipe from durable project storage using HPC scratch.",
    )
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    setup_parser = subparsers.add_parser(
        "setup",
        parents=[_common_parser()],
        help="create editable configuration for a named run",
    )
    setup_parser.add_argument("run_name")
    setup_parser.set_defaults(handler=_setup)

    run_parser = subparsers.add_parser(
        "run",
        parents=[_common_parser(), _scratch_parser()],
        help="run or resume a named workflow",
    )
    run_parser.add_argument("run_name")
    run_parser.add_argument(
        "--restart",
        action="store_true",
        help="discard this run's existing scratch state before starting",
    )
    run_parser.add_argument(
        "--only-fire",
        action="store_true",
        help=(
            "run FIRE and required upstream dependencies while skipping cCRE, "
            "SV, and HLA targets"
        ),
    )
    run_parser.add_argument(
        "--snakemake-executable",
        default=None,
        help=(
            "Snakemake executable (default: LONGPIPE_SNAKEMAKE, the executable "
            "installed alongside Longpipe, or PATH)"
        ),
    )
    run_parser.set_defaults(handler=_run)

    clean_parser = subparsers.add_parser(
        "clean",
        parents=[_common_parser(), _scratch_parser()],
        help="remove retained scratch data for a named run",
    )
    clean_parser.add_argument("run_name")
    clean_parser.set_defaults(handler=_clean)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    try:
        return int(args.handler(args))
    except LongpipeError as error:
        print(f"longpipe: error: {error}", file=sys.stderr)
        return 2
    except OSError as error:
        print(f"longpipe: error: {error}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print("longpipe: interrupted", file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
