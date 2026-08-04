from __future__ import annotations

import contextlib
import io
import json
import os
import stat
import sys
import tempfile
import unittest
from importlib.metadata import version
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "src"))

from longpipe import cli  # noqa: E402


@contextlib.contextmanager
def working_directory(path: Path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


@contextlib.contextmanager
def environment(**values: str):
    previous = {key: os.environ.get(key) for key in values}
    os.environ.update(values)
    try:
        yield
    finally:
        for key, value in previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


class LongpipeCliTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.project = self.root / "project"
        self.scratch = self.root / "scratch"
        self.work_one = self.project / "study-one"
        self.work_two = self.project / "study-two"
        self.work_one.mkdir(parents=True)
        self.work_two.mkdir()
        self.fake_snakemake = self.root / "fake-snakemake"
        self.fake_snakemake.write_text(
            """#!/usr/bin/env python3
import json
import os
import sys
from pathlib import Path

scratch = Path(os.environ["LONGPIPE_SCRATCH_DIR"])
(scratch / "logs" / "rules").mkdir(parents=True, exist_ok=True)
(scratch / "logs" / "rules" / "fake.log").write_text("diagnostic\\n")
(scratch / "benchmarks").mkdir(parents=True, exist_ok=True)
(scratch / "benchmarks" / "fake.tsv").write_text("s\\tmax_rss\\n1\\t1\\n")
(scratch / ".snakemake" / "log").mkdir(parents=True, exist_ok=True)
(scratch / ".snakemake" / "log" / "workflow.log").write_text("metadata\\n")

counter_path = scratch / "resume-counter"
count = int(counter_path.read_text()) + 1 if counter_path.exists() else 1
counter_path.write_text(str(count))

results = scratch / "results"
(results / "fire").mkdir(parents=True, exist_ok=True)
(results / "fire" / "result.txt").write_text(f"invocation {count}\\n")
wgs = results / "wgs" / "sample"
wgs.mkdir(parents=True, exist_ok=True)
(wgs / "sample.bam").write_text("bam\\n")
(wgs / "outputs.json").write_text(json.dumps({
    "haplotagged_bam": str(wgs / "sample.bam"),
    "sample_name": "sample",
}))
if os.environ.get("FAKE_LONGPIPE_RESULT") == "failure":
    print("simulated workflow failure")
    sys.exit(7)
print("simulated workflow success")
""",
            encoding="utf-8",
        )
        self.fake_snakemake.chmod(
            self.fake_snakemake.stat().st_mode | stat.S_IXUSR
        )
        self.base_environment = {
            "LONGPIPE_WORKFLOW_ROOT": str(REPOSITORY_ROOT),
        }

    def tearDown(self):
        self.temporary.cleanup()

    def test_version_uses_installed_package_metadata(self):
        stdout = io.StringIO()
        with (
            contextlib.redirect_stdout(stdout),
            self.assertRaisesRegex(SystemExit, "0"),
        ):
            cli.main(["--version"])

        self.assertEqual(stdout.getvalue().strip(), version("longpipe"))

    def invoke(self, arguments: list[str]) -> tuple[int, str, str]:
        stdout = io.StringIO()
        stderr = io.StringIO()
        with (
            environment(**self.base_environment),
            contextlib.redirect_stdout(stdout),
            contextlib.redirect_stderr(stderr),
        ):
            result = cli.main(arguments)
        return result, stdout.getvalue(), stderr.getvalue()

    def setup_run(self, working: Path, run_name: str = "test-run") -> Path:
        with working_directory(working):
            result, _, error = self.invoke(
                [
                    "setup",
                    run_name,
                    "--project-root",
                    str(self.project),
                ]
            )
        self.assertEqual(result, 0, error)
        return working / run_name

    def run_arguments(self, run_name: str = "test-run") -> list[str]:
        return [
            "run",
            run_name,
            "--project-root",
            str(self.project),
            "--scratch-root",
            str(self.scratch),
            "--snakemake-executable",
            str(self.fake_snakemake),
        ]

    def test_setup_creates_templates_and_refuses_overwrite(self):
        run_directory = self.setup_run(self.work_one)
        config = (run_directory / "config.yaml").read_text()
        self.assertIn("runsheet: runsheet.csv", config)
        self.assertNotIn("\nname:", config)
        self.assertNotIn("\noutdir:", config)
        self.assertTrue((run_directory / "runsheet.csv").is_file())

        with working_directory(self.work_one):
            result, _, error = self.invoke(
                ["setup", "test-run", "--project-root", str(self.project)]
            )
        self.assertEqual(result, 2)
        self.assertIn("already exists", error)

    def test_setup_rejects_work_outside_project_root(self):
        outside = self.root / "outside"
        outside.mkdir()
        with working_directory(outside):
            result, _, error = self.invoke(
                ["setup", "test-run", "--project-root", str(self.project)]
            )
        self.assertEqual(result, 2)
        self.assertIn("must be within project root", error)

    def test_success_snapshots_resumes_restarts_and_cleans(self):
        run_directory = self.setup_run(self.work_one)
        arguments = self.run_arguments()

        with working_directory(self.work_one):
            result, output, error = self.invoke(arguments)
        self.assertEqual(result, 0, error)
        self.assertIn("Run completed", output)
        self.assertEqual(
            (run_directory / "results" / "fire" / "result.txt").read_text(),
            "invocation 1\n",
        )
        published_wgs = run_directory / "results" / "wgs" / "sample"
        published_manifest = json.loads(
            (published_wgs / "outputs.json").read_text()
        )
        self.assertEqual(
            published_manifest["haplotagged_bam"],
            str(published_wgs / "sample.bam"),
        )
        self.assertEqual(published_manifest["sample_name"], "sample")
        scratch_manifest = json.loads(
            (
                self.scratch
                / "test-run"
                / "results"
                / "wgs"
                / "sample"
                / "outputs.json"
            ).read_text()
        )
        self.assertEqual(
            scratch_manifest["haplotagged_bam"],
            str(
                self.scratch
                / "test-run"
                / "results"
                / "wgs"
                / "sample"
                / "sample.bam"
            ),
        )
        self.assertTrue((run_directory / "logs" / "rules" / "fake.log").is_file())
        self.assertTrue((run_directory / "benchmarks" / "fake.tsv").is_file())
        self.assertTrue((self.scratch / "test-run").is_dir())

        with working_directory(self.work_one):
            result, _, error = self.invoke(arguments)
        self.assertEqual(result, 0, error)
        self.assertEqual(
            (run_directory / "results" / "fire" / "result.txt").read_text(),
            "invocation 2\n",
        )

        invocation_directories = sorted(
            (run_directory / "run-metadata" / "invocations").iterdir()
        )
        self.assertEqual(len(invocation_directories), 2)
        for invocation in invocation_directories:
            self.assertTrue((invocation / "config.yaml").is_file())
            self.assertTrue((invocation / "runsheet.csv").is_file())
            self.assertTrue((invocation / "status.json").is_file())

        with working_directory(self.work_one):
            result, _, error = self.invoke(arguments + ["--restart"])
        self.assertEqual(result, 0, error)
        self.assertEqual(
            (run_directory / "results" / "fire" / "result.txt").read_text(),
            "invocation 1\n",
        )

        with working_directory(self.work_one):
            result, _, error = self.invoke(
                [
                    "clean",
                    "test-run",
                    "--project-root",
                    str(self.project),
                    "--scratch-root",
                    str(self.scratch),
                ]
            )
        self.assertEqual(result, 0, error)
        self.assertFalse((self.scratch / "test-run").exists())
        self.assertTrue((run_directory / "results").is_dir())

    def test_failure_publishes_diagnostics_but_not_results(self):
        run_directory = self.setup_run(self.work_one, "failed-run")
        with working_directory(self.work_one), environment(
            FAKE_LONGPIPE_RESULT="failure"
        ):
            result, _, error = self.invoke(self.run_arguments("failed-run"))

        self.assertEqual(result, 7)
        self.assertIn("Run failed with exit code 7", error)
        self.assertFalse((run_directory / "results").exists())
        self.assertTrue((run_directory / "logs" / "rules" / "fake.log").is_file())
        status = json.loads(
            (run_directory / "run-metadata" / "status.json").read_text()
        )
        self.assertEqual(status["status"], "failed")
        self.assertEqual(status["return_code"], 7)
        self.assertTrue((self.scratch / "failed-run" / "results").is_dir())

    def test_scratch_ownership_prevents_cross_project_collision(self):
        self.setup_run(self.work_one, "shared-name")
        self.setup_run(self.work_two, "shared-name")
        with working_directory(self.work_one):
            result, _, error = self.invoke(self.run_arguments("shared-name"))
        self.assertEqual(result, 0, error)

        with working_directory(self.work_two):
            result, _, error = self.invoke(self.run_arguments("shared-name"))
        self.assertEqual(result, 2)
        self.assertIn("belongs to", error)

    def test_active_lock_is_rejected(self):
        self.setup_run(self.work_one, "locked-run")
        run_directory = self.work_one / "locked-run"
        self.scratch.mkdir()
        with cli._run_lock(self.scratch, "locked-run", run_directory):
            with working_directory(self.work_one):
                result, _, error = self.invoke(self.run_arguments("locked-run"))
        self.assertEqual(result, 2)
        self.assertIn("already active", error)

    def test_symlinked_scratch_path_is_rejected(self):
        self.setup_run(self.work_one, "linked-scratch")
        self.scratch.mkdir()
        outside_scratch = self.root / "outside-scratch"
        outside_scratch.mkdir()
        (self.scratch / "linked-scratch").symlink_to(
            outside_scratch, target_is_directory=True
        )

        with working_directory(self.work_one):
            result, _, error = self.invoke(self.run_arguments("linked-scratch"))
        self.assertEqual(result, 2)
        self.assertIn("symlinked scratch path", error)

    def test_run_directory_cannot_escape_project_through_symlink(self):
        outside_run = self.root / "outside-run"
        outside_run.mkdir()
        (self.work_one / "linked-run").symlink_to(
            outside_run, target_is_directory=True
        )

        with working_directory(self.work_one):
            result, _, error = self.invoke(self.run_arguments("linked-run"))
        self.assertEqual(result, 2)
        self.assertIn("resolves outside project root", error)

    def test_cli_managed_config_fields_are_rejected_before_scratch_creation(self):
        run_directory = self.setup_run(self.work_one, "reserved")
        config_path = run_directory / "config.yaml"
        config_path.write_text(f"name: wrong\n{config_path.read_text()}")

        with working_directory(self.work_one):
            result, _, error = self.invoke(self.run_arguments("reserved"))
        self.assertEqual(result, 2)
        self.assertIn("managed by the CLI", error)
        self.assertFalse((self.scratch / "reserved").exists())


if __name__ == "__main__":
    unittest.main()
