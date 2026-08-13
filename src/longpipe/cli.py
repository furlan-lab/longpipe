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
SETUP_RUNSHEET_COLUMNS = (
    "core_name",
    "lab_name",
    "run",
    "cell",
    "barcode",
    "hifi_bam",
    "input_stage",
    "pbwgs_outputs_json",
    "sex",
    "hia5",
    "fiberseq_input",
    "notes",
)


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


def _parse_runsheet_rows(value: str | None) -> tuple[int, ...] | None:
    if value is None:
        return None
    if not value.strip():
        raise LongpipeError("--rows must not be empty")

    selected: set[int] = set()
    for component in value.split(","):
        token = component.strip()
        if not token:
            raise LongpipeError(
                f"invalid --rows value {value!r}: empty selection component"
            )
        bounds = token.split(":")
        if len(bounds) > 2 or any(
            not re.fullmatch(r"[0-9]+", bound.strip()) for bound in bounds
        ):
            raise LongpipeError(
                f"invalid --rows component {token!r}: expected a row number "
                "or inclusive START:END range"
            )
        start = int(bounds[0])
        end = int(bounds[-1])
        if start < 1 or end < 1:
            raise LongpipeError(
                f"invalid --rows component {token!r}: row numbers are 1-based"
            )
        if start > end:
            raise LongpipeError(
                f"invalid --rows range {token!r}: start exceeds end"
            )
        selected.update(range(start, end + 1))
    return tuple(sorted(selected))


def _selected_runsheet_records(
    runsheet_path: Path, selected_rows: Sequence[int]
) -> list[list[str]]:
    try:
        with runsheet_path.open(newline="", encoding="utf-8") as handle:
            records = list(csv.reader(handle))
    except (OSError, csv.Error) as error:
        raise LongpipeError(
            f"could not read runsheet {runsheet_path}: {error}"
        ) from error
    if not records:
        raise LongpipeError(f"runsheet is empty: {runsheet_path}")

    header = records[0]
    data_rows = [record for record in records[1:] if any(record)]
    row_count = len(data_rows)
    unavailable = [row for row in selected_rows if row > row_count]
    if unavailable:
        raise LongpipeError(
            f"--rows selects row {unavailable[0]}, but {runsheet_path} has only "
            f"{row_count} data row{'s' if row_count != 1 else ''}"
        )
    return [header, *(data_rows[row - 1] for row in selected_rows)]


def _snapshot_inputs(
    run_directory: Path,
    scratch_directory: Path,
    invocation_id: str,
    selected_runsheet: Sequence[Sequence[str]] | None = None,
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
    if selected_runsheet is None:
        shutil.copy2(runsheet_path, runsheet_snapshot)
    else:
        with runsheet_snapshot.open("w", newline="", encoding="utf-8") as handle:
            csv.writer(handle, lineterminator="\n").writerows(selected_runsheet)

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


def _replace_selected_sample_trees(
    durable_samples: Path,
    sample_names: Sequence[str],
    build: Callable[[Path], None],
) -> None:
    durable_samples.parent.mkdir(parents=True, exist_ok=True)
    stage = durable_samples.parent / (
        f".{durable_samples.name}.rows-stage-{uuid.uuid4().hex}"
    )
    backups = durable_samples.parent / (
        f".{durable_samples.name}.rows-backup-{uuid.uuid4().hex}"
    )
    published: list[str] = []
    pending: tuple[str, bool] | None = None
    publication_complete = False
    rollback_complete = False

    def remove(path: Path) -> None:
        if path.is_symlink() or path.is_file():
            path.unlink()
        elif path.is_dir():
            shutil.rmtree(path)

    try:
        build(stage)
        durable_samples.mkdir(exist_ok=True)
        backups.mkdir()
        for sample in sample_names:
            source = stage / sample
            destination = durable_samples / sample
            backup = backups / sample
            had_destination = destination.exists() or destination.is_symlink()
            pending = (sample, had_destination)
            if had_destination:
                os.replace(destination, backup)
            os.replace(source, destination)
            published.append(sample)
            pending = None
        publication_complete = True
    except BaseException:
        if pending is not None:
            sample, had_destination = pending
            source = stage / sample
            destination = durable_samples / sample
            backup = backups / sample
            if had_destination:
                if backup.exists() or backup.is_symlink():
                    remove(destination)
                    os.replace(backup, destination)
            elif not source.exists():
                remove(destination)
        for sample in reversed(published):
            destination = durable_samples / sample
            backup = backups / sample
            remove(destination)
            if backup.exists() or backup.is_symlink():
                os.replace(backup, destination)
        rollback_complete = True
        raise
    finally:
        remove(stage)
        if publication_complete or rollback_complete:
            remove(backups)


def _publish_sample_results(
    scratch_results: Path,
    run_directory: Path,
    runsheet_snapshot: Path,
    config: dict,
    preserve_unselected: bool = False,
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

    def build(stage: Path) -> None:
        _build_sample_results(
            stage,
            scratch_results,
            durable_samples,
            sample_names,
            set(ccre_peaks),
            fire_version,
        )
    if preserve_unselected:
        _replace_selected_sample_trees(
            durable_samples,
            sample_names,
            build,
        )
    else:
        _replace_built_tree(
            durable_samples,
            build,
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


def _finalize_per_sample_publication(
    run_directory: Path,
    runsheet_snapshot: Path,
    preserve_unselected: bool,
) -> None:
    """Remove obsolete legacy results after every selected sample published."""

    if not preserve_unselected:
        selected = set(_sample_labels(runsheet_snapshot))
        durable_samples = run_directory / "samples"
        if durable_samples.is_dir():
            for entry in durable_samples.iterdir():
                if entry.name not in selected:
                    if entry.is_dir() and not entry.is_symlink():
                        shutil.rmtree(entry)
                    else:
                        entry.unlink()

    legacy_results = run_directory / "results"
    if legacy_results.is_symlink() or legacy_results.is_file():
        legacy_results.unlink()
    elif legacy_results.is_dir():
        shutil.rmtree(legacy_results)


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


def _publish_fiberseq_metadata(
    scratch_directory: Path,
    run_directory: Path,
    durable_invocation: Path,
    runsheet_snapshot: Path,
    config: dict,
    preserve_unselected: bool = False,
) -> None:
    source = (
        scratch_directory
        / "intermediate"
        / "fiberseq"
        / "input-detection"
    )
    selected_samples = set(_sample_labels(runsheet_snapshot))
    sample_names = selected_samples
    if config.get("hia5_gates_fiberseq", True):
        with runsheet_snapshot.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        sample_names = {
            row.get("lab_name", "").strip()
            or row.get("core_name", "").strip()
            for row in rows
            if row.get("hia5", "").strip().lower()
            in {"yes", "y", "true", "t", "1"}
        }

    def replace(destination: Path, preserve: bool) -> None:
        def build(stage: Path) -> None:
            if preserve and destination.is_dir():
                shutil.copytree(destination, stage)
                for sample in selected_samples:
                    (stage / f"{sample}.json").unlink(missing_ok=True)
            else:
                stage.mkdir()
            for sample in sample_names:
                report = source / f"{sample}.json"
                if report.is_file():
                    shutil.copy2(report, stage / report.name)

        _replace_built_tree(destination, build)

    replace(
        run_directory / "run-metadata" / "fiberseq-input",
        preserve_unselected,
    )
    replace(durable_invocation / "fiberseq-input", False)


def _publish_reanalysis_metadata(
    scratch_directory: Path,
    run_directory: Path,
    durable_invocation: Path,
    runsheet_snapshot: Path,
    preserve_unselected: bool = False,
) -> None:
    source = (
        scratch_directory
        / "intermediate"
        / "reanalysis"
        / "input-validation"
    )
    with runsheet_snapshot.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    selected_samples = set(_sample_labels(runsheet_snapshot))
    sample_names = {
        row.get("lab_name", "").strip()
        or row.get("core_name", "").strip()
        for row in rows
        if row.get("input_stage", "").strip().lower()
        in {"aligned", "pbwgs"}
    }

    def replace(destination: Path, preserve: bool) -> None:
        def build(stage: Path) -> None:
            if preserve and destination.is_dir():
                shutil.copytree(destination, stage)
                for sample in selected_samples:
                    (stage / f"{sample}.json").unlink(missing_ok=True)
            else:
                stage.mkdir()
            for sample in sample_names:
                report = source / f"{sample}.json"
                if report.is_file():
                    shutil.copy2(report, stage / report.name)

        _replace_built_tree(destination, build)

    replace(
        run_directory / "run-metadata" / "reanalysis-input",
        preserve_unselected,
    )
    replace(durable_invocation / "reanalysis-input", False)


def _publish_diagnostics(
    scratch_directory: Path,
    run_directory: Path,
    durable_invocation: Path,
    runsheet_snapshot: Path,
    config: dict,
    preserve_unselected: bool = False,
) -> None:
    _replace_tree(scratch_directory / "logs", run_directory / "logs")
    _replace_tree(scratch_directory / "benchmarks", run_directory / "benchmarks")
    _publish_snakemake_metadata(scratch_directory, durable_invocation)
    _publish_fiberseq_metadata(
        scratch_directory,
        run_directory,
        durable_invocation,
        runsheet_snapshot,
        config,
        preserve_unselected,
    )
    _publish_reanalysis_metadata(
        scratch_directory,
        run_directory,
        durable_invocation,
        runsheet_snapshot,
        preserve_unselected,
    )


def _snakemake_command(
    executable: str,
    workflow_root: Path,
    config_snapshot: Path,
    targets: Sequence[str] = (),
    dry_run: bool = False,
) -> list[str]:
    command = [
        executable,
        *targets,
        "--snakefile",
        str(workflow_root / "Snakefile"),
        "--profile",
        str(workflow_root / "profiles" / "slurm"),
        "--configfile",
        str(config_snapshot),
    ]
    if dry_run:
        command.append("--dry-run")
    return command


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


def _ft_predict_backend(config: dict) -> str:
    configured_backend = config.get("ft_predict_backend")
    configured_executable = config.get("ft_predict_executable", "ft")
    if configured_backend is None:
        return "external" if configured_executable != "ft" else "conda"
    if configured_backend not in {"conda", "external"}:
        raise LongpipeError(
            "ft_predict_backend must be either 'conda' or 'external'"
        )
    return configured_backend


def _run_ft_preflight(
    executable: Path,
    arguments: Sequence[str],
    environment: dict[str, str],
    run_directory: Path,
) -> str:
    command = [str(executable), *arguments]
    try:
        completed = subprocess.run(
            command,
            cwd=run_directory,
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            errors="replace",
            timeout=15,
            check=False,
        )
    except subprocess.TimeoutExpired as error:
        raise LongpipeError(
            f"external ft preflight timed out: {' '.join(command)}"
        ) from error
    except OSError as error:
        raise LongpipeError(
            f"could not run external ft executable {executable}: {error}"
        ) from error

    output = completed.stdout.strip()
    if completed.returncode != 0:
        detail = f": {output}" if output else ""
        raise LongpipeError(
            f"external ft preflight failed with exit code "
            f"{completed.returncode} while running {' '.join(command)}{detail}"
        )
    return output


def _validate_ft_predict_config(
    config: dict, run_directory: Path
) -> dict[str, str]:
    backend = _ft_predict_backend(config)
    if backend == "conda":
        return {"backend": backend}

    configured_executable = config.get("ft_predict_executable")
    if (
        not isinstance(configured_executable, str)
        or not configured_executable.strip()
    ):
        raise LongpipeError(
            "ft_predict_executable must be set when ft_predict_backend is external"
        )
    expanded = os.path.expandvars(
        os.path.expanduser(configured_executable.strip())
    )
    if os.sep not in expanded:
        raise LongpipeError(
            "external ft_predict_executable must be an explicit absolute or "
            "run-relative path, not a PATH-only command"
        )
    executable = _path(expanded, run_directory)
    if not executable.is_file():
        raise LongpipeError(f"external ft executable does not exist: {executable}")
    if not os.access(executable, os.X_OK):
        raise LongpipeError(f"external ft executable is not executable: {executable}")

    environment = os.environ.copy()
    for variable in (
        "LIBTORCH",
        "LIBTORCH_CXX11_ABI",
        "LD_LIBRARY_PATH",
        "DYLD_LIBRARY_PATH",
        "LD_PRELOAD",
    ):
        environment.pop(variable, None)

    version_output = _run_ft_preflight(
        executable,
        ("--version",),
        environment,
        run_directory,
    )
    version_match = re.search(
        r"\bfibertools-rs\s+v?(\d+\.\d+\.\d+)\b", version_output
    )
    if version_match is None:
        raise LongpipeError(
            f"could not determine fibertools-rs version from {executable}: "
            f"{version_output!r}"
        )
    version = version_match.group(1)

    _run_ft_preflight(
        executable,
        ("predict-m6a", "--help"),
        environment,
        run_directory,
    )

    expected_version = config.get("ft_predict_expected_version")
    result = {
        "backend": backend,
        "executable": str(executable),
        "version": version,
    }
    if expected_version is not None:
        if (
            not isinstance(expected_version, str)
            or not expected_version.strip()
        ):
            raise LongpipeError(
                "ft_predict_expected_version must be a non-empty version string"
            )
        normalized_expected = expected_version.strip().removeprefix("v")
        if normalized_expected != version:
            raise LongpipeError(
                f"external ft version mismatch: expected {normalized_expected}, "
                f"found {version} at {executable}"
            )
        result["expected_version"] = normalized_expected
    return result


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
        with runsheet_template.open(newline="", encoding="utf-8") as handle:
            header = next(csv.reader(handle))
    except (OSError, csv.Error, StopIteration) as error:
        raise LongpipeError(
            f"could not read packaged runsheet template {runsheet_template}: {error}"
        ) from error
    if tuple(header) != SETUP_RUNSHEET_COLUMNS:
        raise LongpipeError(
            "packaged runsheet template is obsolete or incompatible; expected "
            f"{','.join(SETUP_RUNSHEET_COLUMNS)}, found {','.join(header)}"
        )
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
    if args.dry_run and args.restart:
        raise LongpipeError("--dry-run cannot be combined with --restart")

    run_directory = _run_directory(args.project_root, args.run_name)
    if not run_directory.is_dir():
        raise LongpipeError(f"run directory does not exist: {run_directory}")
    run_config, runsheet_path = _load_run_config(
        run_directory / "config.yaml", run_directory
    )
    row_selection = _parse_runsheet_rows(args.rows)
    selected_runsheet = (
        _selected_runsheet_records(runsheet_path, row_selection)
        if row_selection is not None
        else None
    )

    scratch_root, scratch_directory = _scratch_directory(
        args.scratch_root, args.run_name
    )
    workflow_root = _workflow_root()
    snakemake_executable = _resolve_snakemake_executable(
        args.snakemake_executable
    )
    ft_predict_runtime = _validate_ft_predict_config(run_config, run_directory)
    if ft_predict_runtime["backend"] == "external":
        print(
            f"Validated external ft {ft_predict_runtime['version']}: "
            f"{ft_predict_runtime['executable']}"
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
            run_directory,
            scratch_directory,
            invocation_id,
            selected_runsheet,
        )
        config, _ = _load_run_config(config_snapshot, config_snapshot.parent)
        run_mode = "fire_only" if args.only_fire else "full"
        targets = ("fire_only",) if args.only_fire else ("all",)
        command = _snakemake_command(
            snakemake_executable,
            workflow_root,
            config_snapshot,
            targets,
            dry_run=args.dry_run,
        )
        status: dict[str, object] = {
            "command": command,
            "dry_run": args.dry_run,
            "hostname": socket.gethostname(),
            "invocation_id": invocation_id,
            "pid": os.getpid(),
            "run_directory": str(run_directory),
            "mode": run_mode,
            "run_name": args.run_name,
            "scratch_directory": str(scratch_directory),
            "started_at": _now(),
            "status": "dry_run_running" if args.dry_run else "running",
            "ft_predict": ft_predict_runtime,
        }
        if row_selection is not None:
            status["runsheet_row_selector"] = args.rows
            status["runsheet_rows"] = list(row_selection)
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
        if row_selection is not None:
            mode_label += f" for runsheet rows {args.rows}"
        action = "Dry-running" if args.dry_run else "Running"
        print(f"{action} {mode_label} for {args.run_name!r} in {scratch_directory}")
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
        if args.dry_run:
            status["status"] = (
                "dry_run_success" if return_code == 0 else "dry_run_failed"
            )
            _record_status(
                status,
                scratch_directory,
                run_directory,
                durable_invocation,
            )
            if return_code == 0:
                print(
                    "Dry-run completed; no jobs were submitted and no results "
                    "were published."
                )
            else:
                print(
                    f"Dry-run failed with exit code {return_code}; no results "
                    "were published.",
                    file=sys.stderr,
                )
            return return_code

        status["status"] = "finalizing" if return_code == 0 else "failed"
        _record_status(status, scratch_directory, run_directory, durable_invocation)

        try:
            _publish_diagnostics(
                scratch_directory,
                run_directory,
                durable_invocation,
                runsheet_snapshot,
                config,
                row_selection is not None,
            )
            if return_code == 0:
                if not args.only_fire:
                    _finalize_per_sample_publication(
                        run_directory,
                        runsheet_snapshot,
                        preserve_unselected=row_selection is not None,
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
        "-n",
        "--dry-run",
        action="store_true",
        help=(
            "preview the selected workflow without submitting jobs or publishing "
            "results"
        ),
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
        "--rows",
        metavar="SELECTION",
        help=(
            "run only selected 1-based runsheet data rows; accepts comma-separated "
            "rows and inclusive ranges such as 1,3,4:8"
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
