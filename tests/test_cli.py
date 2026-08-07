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
from unittest import mock


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
import csv
import json
import os
import sys
from pathlib import Path

scratch = Path(os.environ["LONGPIPE_SCRATCH_DIR"])
fire_only = "fire_only" in sys.argv[1:]
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
with Path(os.environ["LONGPIPE_RUNSHEET"]).open(newline="") as handle:
    rows = list(csv.DictReader(handle))
samples = [(row["lab_name"].strip() or row["core_name"].strip()) for row in rows]
fire_samples = [
    sample
    for sample, row in zip(samples, rows)
    if row["hia5"].strip().lower() in {"yes", "y", "true", "t", "1"}
]

fire = results / "fire"
fire.mkdir(parents=True, exist_ok=True)
(fire / "eligibility.tsv").write_text(
    "sample\\tmean_depth\\tminimum_mean_depth\\teligible\\treason\\n"
    + "".join(
        f"{sample}\\t20\\t8\\tyes\\teligible\\n" for sample in fire_samples
    )
)
for sample in fire_samples:
    pileups = fire / "pileups"
    pileups.mkdir(parents=True, exist_ok=True)
    (pileups / f"{sample}-fire-v0.1.2-pileup.bed.gz").write_text(
        f"{'fire-only ' if fire_only else ''}invocation {count}\\n"
    )
    (pileups / f"{sample}-fire-v0.1.2-pileup.bed.gz.tbi").write_text(
        "index\\n"
    )
    bigwigs = fire / "bigwigs"
    bigwigs.mkdir(parents=True, exist_ok=True)
    (bigwigs / f"{sample}-fire-v0.1.2-pileup_shrunk.bw").write_text(
        "bigwig\\n"
    )
    tracks = fire / "haplotypes" / "tracks"
    tracks.mkdir(parents=True, exist_ok=True)
    for track_name in (
        f"{sample}-fire-v0.1.2-H1-percent-accessible.bw",
        f"{sample}-fire-v0.1.2-H2-percent-accessible.bw",
        f"{sample}-fire-v0.1.2-H1-minus-H2-percent-accessible.bw",
    ):
        (tracks / track_name).write_text("haplotype track\\n")
    differences = fire / "haplotypes" / "differences"
    differences.mkdir(parents=True, exist_ok=True)
    differences_name = f"{sample}-fire-v0.1.2-haplotype-differences.bed.gz"
    (differences / differences_name).write_text("haplotype differences\\n")
    (differences / f"{differences_name}.tbi").write_text("index\\n")
    (
        differences
        / f"{sample}-fire-v0.1.2-haplotype-differences-summary.tsv"
    ).write_text("metric\\tvalue\\n")
    if not fire_only:
        for peak_set in ("brain_atac", "marrow_atac"):
            ccre = fire / "ccre" / peak_set
            ccre.mkdir(parents=True, exist_ok=True)
            for statistic in ("max", "median", "sum"):
                (ccre / f"{sample}_scATAC_{statistic}_scores.bed").write_text(
                    "ccre\\n"
                )

if not fire_only:
    for sample in samples:
        wgs = results / "wgs" / sample
        wgs.mkdir(parents=True, exist_ok=True)
        (wgs / f"{sample}.bam").write_text("bam\\n")
        (wgs / "outputs.json").write_text(json.dumps({
            "haplotagged_bam": str(wgs / f"{sample}.bam"),
            "nested": {"files": [str(wgs / f"{sample}.bam")]},
            "sample_name": sample,
        }))

        sv = results / "sv"
        sv.mkdir(parents=True, exist_ok=True)
        (sv / f"{sample}.sawfish.vcf.gz").write_text("vcf\\n")
        (sv / f"{sample}.sawfish.vcf.gz.tbi").write_text("index\\n")
        (sv / f"{sample}.supporting_reads.json.gz").write_text("reads\\n")
        svtopo = sv / f"{sample}.svtopo"
        svtopo.mkdir(parents=True, exist_ok=True)
        (svtopo / "summary.tsv").write_text("svtopo\\n")

        hla = results / "hla"
        hla.mkdir(parents=True, exist_ok=True)
        (hla / f"{sample}_hifihla_report.json").write_text("{}\\n")
        (hla / f"{sample}_hifihla_report.tsv").write_text("hla\\n")
        (hla / f"{sample}_hifihla_summary.tsv").write_text("hla\\n")
if os.environ.get("FAKE_LONGPIPE_UNEXPECTED") == "yes":
    (results / "unrouted").mkdir(parents=True, exist_ok=True)
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

    def test_run_uses_snakemake_installed_beside_longpipe(self):
        run_directory = self.setup_run(self.work_one, "sibling-snakemake")
        sibling_directory = self.root / "pipx" / "bin"
        sibling_directory.mkdir(parents=True)
        sibling_snakemake = sibling_directory / "snakemake"
        sibling_snakemake.write_bytes(self.fake_snakemake.read_bytes())
        sibling_snakemake.chmod(self.fake_snakemake.stat().st_mode)

        arguments = [
            "run",
            "sibling-snakemake",
            "--project-root",
            str(self.project),
            "--scratch-root",
            str(self.scratch),
        ]
        with (
            mock.patch.object(
                cli.sys, "executable", str(sibling_directory / "python")
            ),
            environment(LONGPIPE_SNAKEMAKE=""),
            working_directory(self.work_one),
        ):
            result, _, error = self.invoke(arguments)

        self.assertEqual(result, 0, error)
        status = json.loads(
            (run_directory / "run-metadata" / "status.json").read_text()
        )
        self.assertEqual(status["command"][0], str(sibling_snakemake))

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
        sample = "Brain_1"
        published_sample = run_directory / "samples" / sample
        published_pileup = (
            published_sample
            / "fire"
            / "pileups"
            / f"{sample}-fire-v0.1.2-pileup.bed.gz"
        )
        self.assertEqual(
            published_pileup.read_text(),
            "invocation 1\n",
        )
        published_wgs = published_sample / "pbWGS"
        published_manifest = json.loads(
            (published_wgs / "outputs.json").read_text()
        )
        self.assertEqual(
            Path(published_manifest["haplotagged_bam"]).resolve(),
            (published_wgs / f"{sample}.bam").resolve(),
        )
        self.assertEqual(
            [Path(path).resolve() for path in published_manifest["nested"]["files"]],
            [(published_wgs / f"{sample}.bam").resolve()],
        )
        self.assertEqual(published_manifest["sample_name"], sample)
        scratch_manifest = json.loads(
            (
                self.scratch
                / "test-run"
                / "results"
                / "wgs"
                / sample
                / "outputs.json"
            ).read_text()
        )
        self.assertEqual(
            Path(scratch_manifest["haplotagged_bam"]).resolve(),
            (
                self.scratch
                / "test-run"
                / "results"
                / "wgs"
                / sample
                / f"{sample}.bam"
            ).resolve(),
        )
        self.assertTrue(
            (published_sample / "sawfish" / f"{sample}.sawfish.vcf.gz").is_file()
        )
        self.assertTrue((published_sample / "svtopo" / "summary.tsv").is_file())
        self.assertTrue(
            (published_sample / "hla" / f"{sample}_hifihla_report.json").is_file()
        )
        self.assertTrue(
            (
                published_sample
                / "fire"
                / "ccre"
                / "brain_atac"
                / f"{sample}_scATAC_max_scores.bed"
            ).is_file()
        )
        self.assertTrue(
            (
                published_sample
                / "fire"
                / "haplotypes"
                / "tracks"
                / f"{sample}-fire-v0.1.2-H1-percent-accessible.bw"
            ).is_file()
        )
        self.assertTrue(
            (
                published_sample
                / "fire"
                / "haplotypes"
                / "differences"
                / f"{sample}-fire-v0.1.2-haplotype-differences.bed.gz"
            ).is_file()
        )
        self.assertFalse((run_directory / "samples" / "BT54" / "fire").exists())
        self.assertTrue(
            (run_directory / "run-metadata" / "fire-eligibility.tsv").is_file()
        )
        self.assertFalse((run_directory / "results").exists())
        self.assertTrue((run_directory / "logs" / "rules" / "fake.log").is_file())
        self.assertTrue((run_directory / "benchmarks" / "fake.tsv").is_file())
        self.assertTrue((self.scratch / "test-run").is_dir())

        legacy_results = run_directory / "results"
        legacy_results.mkdir()
        (legacy_results / "old-result.txt").write_text("old\n")
        with working_directory(self.work_one):
            result, _, error = self.invoke(arguments)
        self.assertEqual(result, 0, error)
        self.assertEqual(published_pileup.read_text(), "invocation 2\n")
        self.assertFalse(legacy_results.exists())

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
        self.assertEqual(published_pileup.read_text(), "invocation 1\n")

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
        self.assertTrue((run_directory / "samples").is_dir())

    def test_fire_only_updates_fire_and_preserves_other_results(self):
        run_directory = self.setup_run(self.work_one, "fire-only-run")
        arguments = self.run_arguments("fire-only-run")
        with working_directory(self.work_one):
            result, _, error = self.invoke(arguments)
        self.assertEqual(result, 0, error)

        sample = "Brain_1"
        published_sample = run_directory / "samples" / sample
        pileup = (
            published_sample
            / "fire"
            / "pileups"
            / f"{sample}-fire-v0.1.2-pileup.bed.gz"
        )
        ccre = (
            published_sample
            / "fire"
            / "ccre"
            / "brain_atac"
            / f"{sample}_scATAC_max_scores.bed"
        )
        ccre.write_text("preserve this cCRE\n")
        wgs_manifest = published_sample / "pbWGS" / "outputs.json"
        sv_vcf = published_sample / "sawfish" / f"{sample}.sawfish.vcf.gz"
        hla_report = published_sample / "hla" / f"{sample}_hifihla_report.json"
        preserved_outputs = {
            path: path.read_bytes()
            for path in (wgs_manifest, sv_vcf, hla_report)
        }

        with working_directory(self.work_one):
            result, output, error = self.invoke(
                arguments + ["--only-fire", "--restart"]
            )

        self.assertEqual(result, 0, error)
        self.assertIn("FIRE-only run completed", output)
        self.assertEqual(pileup.read_text(), "fire-only invocation 1\n")
        self.assertEqual(ccre.read_text(), "preserve this cCRE\n")
        for path, expected in preserved_outputs.items():
            with self.subTest(path=path):
                self.assertEqual(path.read_bytes(), expected)
        self.assertFalse(
            (self.scratch / "fire-only-run" / "results" / "sv").exists()
        )
        self.assertFalse(
            (self.scratch / "fire-only-run" / "results" / "hla").exists()
        )
        status = json.loads(
            (run_directory / "run-metadata" / "status.json").read_text()
        )
        self.assertEqual(status["mode"], "fire_only")
        self.assertEqual(status["command"][-1], "fire_only")
        self.assertEqual(status["status"], "success")

    def test_failure_publishes_diagnostics_but_not_results(self):
        run_directory = self.setup_run(self.work_one, "failed-run")
        with working_directory(self.work_one), environment(
            FAKE_LONGPIPE_RESULT="failure"
        ):
            result, _, error = self.invoke(self.run_arguments("failed-run"))

        self.assertEqual(result, 7)
        self.assertIn("Run failed with exit code 7", error)
        self.assertFalse((run_directory / "samples").exists())
        self.assertTrue((run_directory / "logs" / "rules" / "fake.log").is_file())
        status = json.loads(
            (run_directory / "run-metadata" / "status.json").read_text()
        )
        self.assertEqual(status["status"], "failed")
        self.assertEqual(status["return_code"], 7)
        self.assertTrue((self.scratch / "failed-run" / "results").is_dir())

    def test_resume_excludes_removed_samples_and_stale_fire(self):
        run_directory = self.setup_run(self.work_one, "changed-run")
        arguments = self.run_arguments("changed-run")
        with working_directory(self.work_one):
            result, _, error = self.invoke(arguments)
        self.assertEqual(result, 0, error)
        self.assertTrue((run_directory / "samples" / "Brain_1" / "fire").is_dir())
        self.assertTrue((run_directory / "samples" / "BT54").is_dir())

        (run_directory / "runsheet.csv").write_text(
            "core_name,lab_name,run,cell,barcode,hifi_bam,sex,hia5\n"
            "BT53,Brain_1,r84265_20260410_231712,1_A01,bc2081,,F,no\n"
        )
        with working_directory(self.work_one):
            result, _, error = self.invoke(arguments)
        self.assertEqual(result, 0, error)
        self.assertFalse((run_directory / "samples" / "Brain_1" / "fire").exists())
        self.assertFalse((run_directory / "samples" / "BT54").exists())
        self.assertTrue(
            (
                self.scratch
                / "changed-run"
                / "results"
                / "fire"
                / "pileups"
                / "Brain_1-fire-v0.1.2-pileup.bed.gz"
            ).is_file()
        )
        self.assertTrue(
            (
                self.scratch
                / "changed-run"
                / "results"
                / "wgs"
                / "BT54"
            ).is_dir()
        )

    def test_publication_failure_preserves_previous_sample_tree(self):
        run_directory = self.setup_run(self.work_one, "publication-run")
        arguments = self.run_arguments("publication-run")
        with working_directory(self.work_one):
            result, _, error = self.invoke(arguments)
        self.assertEqual(result, 0, error)
        pileup = (
            run_directory
            / "samples"
            / "Brain_1"
            / "fire"
            / "pileups"
            / "Brain_1-fire-v0.1.2-pileup.bed.gz"
        )
        self.assertEqual(pileup.read_text(), "invocation 1\n")

        with working_directory(self.work_one), environment(
            FAKE_LONGPIPE_UNEXPECTED="yes"
        ):
            result, _, error = self.invoke(arguments)
        self.assertEqual(result, 2)
        self.assertIn("publication failed", error)
        self.assertEqual(pileup.read_text(), "invocation 1\n")
        status = json.loads(
            (run_directory / "run-metadata" / "status.json").read_text()
        )
        self.assertEqual(status["status"], "publication_failed")

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
