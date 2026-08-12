"""Atomic publication of completed per-sample workflow results."""

from __future__ import annotations

import json
import os
import shutil
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable


class PublicationError(RuntimeError):
    """A completed sample could not be published safely."""


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _remove(path: Path) -> None:
    if path.is_symlink() or path.is_file():
        path.unlink()
    elif path.is_dir():
        shutil.rmtree(path)


def _atomic_write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        temporary.write_text(
            json.dumps(payload, indent=2) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _replace_built_tree(
    destination: Path,
    build: Callable[[Path], None],
) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    stage = destination.parent / f".{destination.name}.stage-{uuid.uuid4().hex}"
    backup = destination.parent / (
        f".{destination.name}.backup-{uuid.uuid4().hex}"
    )
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
        _remove(backup)
    finally:
        _remove(stage)
        if backup.exists() or backup.is_symlink():
            if destination.exists() or destination.is_symlink():
                _remove(backup)
            else:
                os.replace(backup, destination)


def _remove_tree_atomically(destination: Path) -> None:
    if not (destination.exists() or destination.is_symlink()):
        return
    backup = destination.parent / (
        f".{destination.name}.remove-{uuid.uuid4().hex}"
    )
    os.replace(destination, backup)
    try:
        _remove(backup)
    except BaseException:
        if not (destination.exists() or destination.is_symlink()):
            os.replace(backup, destination)
        raise


def _copy_required(source: Path, destination: Path) -> None:
    if not source.exists():
        raise PublicationError(f"expected workflow output is missing: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    if source.is_dir():
        shutil.copytree(source, destination)
    else:
        shutil.copy2(source, destination)


def _rewrite_wgs_manifest(
    manifest: Path,
    scratch_sample: Path,
    durable_sample: Path,
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
        raise PublicationError(
            f"could not read WGS output manifest {manifest}: {error}"
        ) from error
    if not isinstance(payload, dict):
        raise PublicationError(
            f"WGS output manifest is not a JSON object: {manifest}"
        )
    manifest.write_text(
        json.dumps(rewrite(payload), indent=2) + "\n",
        encoding="utf-8",
    )


def read_fire_decision(path: Path, sample: str) -> dict[str, object]:
    try:
        decision = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise PublicationError(
            f"could not read FIRE eligibility decision {path}: {error}"
        ) from error
    if not isinstance(decision, dict) or decision.get("sample") != sample:
        raise PublicationError(
            f"FIRE eligibility decision does not match sample {sample}: {path}"
        )
    if not isinstance(decision.get("eligible"), bool):
        raise PublicationError(
            f"FIRE eligibility decision is not boolean for {sample}: {path}"
        )
    return decision


def _copy_fire_outputs(
    scratch_fire: Path,
    destination: Path,
    sample: str,
    fire_version: str,
    peak_sets: Iterable[str] | None,
) -> None:
    pileup = f"{sample}-fire-{fire_version}-pileup.bed.gz"
    for name in (pileup, f"{pileup}.tbi"):
        _copy_required(
            scratch_fire / "pileups" / name,
            destination / "pileups" / name,
        )

    bigwig = f"{sample}-fire-{fire_version}-pileup_shrunk.bw"
    _copy_required(
        scratch_fire / "bigwigs" / bigwig,
        destination / "bigwigs" / bigwig,
    )

    track_names = (
        f"{sample}-fire-{fire_version}-H1-percent-accessible.bw",
        f"{sample}-fire-{fire_version}-H2-percent-accessible.bw",
        f"{sample}-fire-{fire_version}-H1-minus-H2-percent-accessible.bw",
    )
    for name in track_names:
        _copy_required(
            scratch_fire / "haplotypes" / "tracks" / name,
            destination / "haplotypes" / "tracks" / name,
        )

    differences = f"{sample}-fire-{fire_version}-haplotype-differences.bed.gz"
    difference_names = (
        differences,
        f"{differences}.tbi",
        f"{sample}-fire-{fire_version}-haplotype-differences-summary.tsv",
    )
    for name in difference_names:
        _copy_required(
            scratch_fire / "haplotypes" / "differences" / name,
            destination / "haplotypes" / "differences" / name,
        )

    if peak_sets is None:
        return
    for peak_set in peak_sets:
        for statistic in ("max", "median", "sum"):
            name = f"{sample}_scATAC_{statistic}_scores.bed"
            _copy_required(
                scratch_fire / "ccre" / peak_set / name,
                destination / "ccre" / peak_set / name,
            )


def _publication_status(
    *,
    sample: str,
    mode: str,
    scratch_directory: Path,
    invocation_id: str,
    fire_decision: dict[str, object] | None,
) -> dict[str, object]:
    return {
        "sample": sample,
        "mode": mode,
        "status": "complete",
        "published_at": _now(),
        "scratch_directory": str(scratch_directory),
        "invocation_id": invocation_id,
        "fire_eligibility": fire_decision,
    }


def publish_sample(
    *,
    sample: str,
    scratch_directory: Path,
    run_directory: Path,
    status_path: Path,
    fire_version: str,
    peak_sets: Iterable[str],
    fire_decision_path: Path | None,
    invocation_id: str = "",
) -> None:
    """Publish one completed full-workflow sample atomically."""

    scratch_results = scratch_directory / "results"
    scratch_wgs = scratch_results / "wgs" / sample
    destination = run_directory / "samples" / sample
    decision = (
        read_fire_decision(fire_decision_path, sample)
        if fire_decision_path is not None
        else None
    )

    def build(stage: Path) -> None:
        stage.mkdir()
        _copy_required(scratch_wgs, stage / "pbWGS")
        manifest = stage / "pbWGS" / "outputs.json"
        if not manifest.is_file():
            raise PublicationError(f"WGS output manifest is missing: {scratch_wgs}")
        _rewrite_wgs_manifest(
            manifest,
            scratch_wgs,
            destination / "pbWGS",
        )

        if decision is not None and decision["eligible"]:
            _copy_fire_outputs(
                scratch_results / "fire",
                stage / "fire",
                sample,
                fire_version,
                peak_sets,
            )

        scratch_sv = scratch_results / "sv"
        for name in (
            f"{sample}.sawfish.vcf.gz",
            f"{sample}.sawfish.vcf.gz.tbi",
            f"{sample}.supporting_reads.json.gz",
        ):
            _copy_required(scratch_sv / name, stage / "sawfish" / name)
        _copy_required(
            scratch_sv / f"{sample}.svtopo",
            stage / "svtopo",
        )

        scratch_hla = scratch_results / "hla"
        for suffix in (
            "_hifihla_report.json",
            "_hifihla_report.tsv",
            "_hifihla_summary.tsv",
        ):
            name = f"{sample}{suffix}"
            _copy_required(scratch_hla / name, stage / "hla" / name)

    _replace_built_tree(destination, build)
    _atomic_write_json(
        status_path,
        _publication_status(
            sample=sample,
            mode="full",
            scratch_directory=scratch_directory,
            invocation_id=invocation_id,
            fire_decision=decision,
        ),
    )


def publish_fire_only(
    *,
    sample: str,
    scratch_directory: Path,
    run_directory: Path,
    status_path: Path,
    fire_version: str,
    fire_decision_path: Path,
    invocation_id: str = "",
) -> None:
    """Publish one sample's FIRE core outputs without replacing other results."""

    decision = read_fire_decision(fire_decision_path, sample)
    destination = run_directory / "samples" / sample / "fire"
    if decision["eligible"]:

        def build(stage: Path) -> None:
            stage.mkdir()
            existing_ccre = destination / "ccre"
            if existing_ccre.is_dir():
                shutil.copytree(existing_ccre, stage / "ccre")
            _copy_fire_outputs(
                scratch_directory / "results" / "fire",
                stage,
                sample,
                fire_version,
                None,
            )

        _replace_built_tree(destination, build)
    else:
        _remove_tree_atomically(destination)

    _atomic_write_json(
        status_path,
        _publication_status(
            sample=sample,
            mode="fire_only",
            scratch_directory=scratch_directory,
            invocation_id=invocation_id,
            fire_decision=decision,
        ),
    )
