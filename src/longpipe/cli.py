"""Manage durable Longpipe runs backed by temporary scratch storage."""

from __future__ import annotations

import argparse
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
from typing import Iterator, Sequence

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
    raise LongpipeError(f"could not locate Longpipe workflow resources; checked {checked}")


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
    if working_directory != project_root and project_root not in working_directory.parents:
        raise LongpipeError(
            f"working directory must be within project root {project_root}; "
            f"got {working_directory}"
        )
    return project_root, working_directory


def _run_directory(project_root_value: str, run_name: str) -> Path:
    _validate_run_name(run_name)
    _, working_directory = _validate_working_directory(project_root_value)
    return working_directory / run_name


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
        raise LongpipeError(f"invalid scratch ownership metadata: {owner_path}") from error
    if not isinstance(owner, dict):
        raise LongpipeError(f"invalid scratch ownership metadata: {owner_path}")
    return owner


def _verify_owner(
    scratch_directory: Path, run_name: str, run_directory: Path
) -> dict[str, object]:
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
        raise LongpipeError(f"refusing to remove symlinked scratch path: {scratch_directory}")
    if scratch_directory.parent != scratch_root:
        raise LongpipeError(f"refusing to remove unsafe scratch path: {scratch_directory}")
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


def _replace_tree(source: Path, destination: Path, required: bool = False) -> None:
    if required and not source.is_dir():
        raise LongpipeError(f"expected workflow output directory is missing: {source}")

    destination.parent.mkdir(parents=True, exist_ok=True)
    stage = destination.parent / f".{destination.name}.stage-{uuid.uuid4().hex}"
    backup = destination.parent / f".{destination.name}.backup-{uuid.uuid4().hex}"
    try:
        if source.is_dir():
            shutil.copytree(source, stage)
        else:
            stage.mkdir()

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
    executable: str, workflow_root: Path, config_snapshot: Path
) -> list[str]:
    return [
        executable,
        "--snakefile",
        str(workflow_root / "Snakefile"),
        "--profile",
        str(workflow_root / "profiles" / "slurm"),
        "--configfile",
        str(config_snapshot),
    ]


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
        command = _snakemake_command(
            args.snakemake_executable, workflow_root, config_snapshot
        )
        status: dict[str, object] = {
            "command": command,
            "hostname": socket.gethostname(),
            "invocation_id": invocation_id,
            "pid": os.getpid(),
            "run_directory": str(run_directory),
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

        print(f"Running {args.run_name!r} in {scratch_directory}")
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
                _replace_tree(
                    scratch_directory / "results",
                    run_directory / "results",
                    required=True,
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
            print(f"Run completed; results published to {run_directory / 'results'}")
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
        "--snakemake-executable",
        default=os.environ.get("LONGPIPE_SNAKEMAKE", "snakemake"),
        help="Snakemake executable (default: LONGPIPE_SNAKEMAKE or %(default)s)",
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
    except KeyboardInterrupt:
        print("longpipe: interrupted", file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
