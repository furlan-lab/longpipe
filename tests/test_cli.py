from __future__ import annotations

import contextlib
import csv
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
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(os.environ["LONGPIPE_WORKFLOW_ROOT"]) / "src"))
from longpipe.publication import publish_fire_only, publish_sample

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
reanalysis_validation = (
    scratch / "intermediate" / "reanalysis" / "input-validation"
)
reanalysis_validation.mkdir(parents=True, exist_ok=True)
for sample, row in zip(samples, rows):
    input_stage = row.get("input_stage", "").strip().lower()
    if input_stage in {"aligned", "pbwgs"}:
        (reanalysis_validation / f"{sample}.json").write_text(json.dumps({
            "sample": sample,
            "input_stage": input_stage,
            "validation": "passed",
        }))
fiberseq_detection = (
    scratch / "intermediate" / "fiberseq" / "input-detection"
)
fiberseq_detection.mkdir(parents=True, exist_ok=True)
for sample in fire_samples:
    (fiberseq_detection / f"{sample}.json").write_text(json.dumps({
        "sample": sample,
        "requested_mode": "auto",
        "selected_action": "use_existing",
        "decision_reason": "m6a_and_footprints_detected",
    }))

fire = results / "fire"
fire.mkdir(parents=True, exist_ok=True)
(fire / "eligibility.tsv").write_text(
    "sample\\tmean_depth\\tminimum_mean_depth\\teligible\\treason\\n"
    + "".join(
        f"{sample}\\t20\\t8\\tyes\\teligible\\n" for sample in fire_samples
    )
)
for sample in fire_samples:
    eligibility = scratch / "intermediate" / "fire" / "eligibility"
    eligibility.mkdir(parents=True, exist_ok=True)
    (eligibility / f"{sample}.json").write_text(json.dumps({
        "sample": sample,
        "mean_depth": 20.0,
        "minimum_mean_depth": 8.0,
        "eligible": True,
        "reason": "eligible",
    }))
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
        (wgs / f"{sample}.bam").write_text(f"invocation {count}\\n")
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
if os.environ.get("FAKE_LONGPIPE_RESULT") == "failure":
    print("simulated workflow failure")
    sys.exit(7)

run_directory = Path(os.environ["LONGPIPE_RUN_DIR"])
invocation_id = os.environ.get("LONGPIPE_INVOCATION_ID", "")
if os.environ.get("FAKE_LONGPIPE_UNEXPECTED") == "yes" and len(samples) > 1:
    (
        results / "hla" / f"{samples[1]}_hifihla_report.json"
    ).unlink(missing_ok=True)

if fire_only:
    for sample in fire_samples:
        publish_fire_only(
            sample=sample,
            scratch_directory=scratch,
            run_directory=run_directory,
            status_path=(
                run_directory
                / "run-metadata"
                / "fire-publication"
                / f"{sample}.json"
            ),
            fire_version="v0.1.2",
            fire_decision_path=(
                scratch
                / "intermediate"
                / "fire"
                / "eligibility"
                / f"{sample}.json"
            ),
            invocation_id=invocation_id,
        )
else:
    for sample in samples:
        decision = (
            scratch
            / "intermediate"
            / "fire"
            / "eligibility"
            / f"{sample}.json"
        )
        publish_sample(
            sample=sample,
            scratch_directory=scratch,
            run_directory=run_directory,
            status_path=(
                run_directory
                / "run-metadata"
                / "sample-publication"
                / f"{sample}.json"
            ),
            fire_version="v0.1.2",
            peak_sets=("brain_atac", "marrow_atac"),
            fire_decision_path=decision if sample in fire_samples else None,
            invocation_id=invocation_id,
        )

(run_directory / "run-metadata").mkdir(parents=True, exist_ok=True)
shutil.copy2(
    fire / "eligibility.tsv",
    run_directory / "run-metadata" / "fire-eligibility.tsv",
)
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

    def test_runsheet_row_selector_parses_rows_and_inclusive_ranges(self):
        self.assertEqual(
            cli._parse_runsheet_rows("1,3,4:8"),
            (1, 3, 4, 5, 6, 7, 8),
        )
        self.assertEqual(
            cli._parse_runsheet_rows(" 3, 1, 2:4, 3 "),
            (1, 2, 3, 4),
        )
        self.assertIsNone(cli._parse_runsheet_rows(None))

        invalid = {
            "": "must not be empty",
            "0": "1-based",
            "1,,3": "empty selection component",
            "4:2": "start exceeds end",
            "1:2:3": "expected a row number",
            "one": "expected a row number",
        }
        for value, message in invalid.items():
            with self.subTest(value=value), self.assertRaisesRegex(
                cli.LongpipeError, message
            ):
                cli._parse_runsheet_rows(value)

    def test_ft_predict_backend_preserves_legacy_configuration(self):
        self.assertEqual(
            cli._ft_predict_backend({"ft_predict_executable": "ft"}),
            "conda",
        )
        self.assertEqual(
            cli._ft_predict_backend(
                {"ft_predict_executable": "/site/bin/gpu-ft"}
            ),
            "external",
        )
        self.assertEqual(
            cli._ft_predict_backend(
                {
                    "ft_predict_backend": "conda",
                    "ft_predict_executable": "/site/bin/gpu-ft",
                }
            ),
            "conda",
        )

    def test_external_ft_predict_preflight_uses_clean_runtime(self):
        binary = self.root / "bin" / "ft"
        binary.parent.mkdir()
        binary.write_text(
            """#!/bin/sh
if [ "${LD_LIBRARY_PATH+x}" = x ] || [ "${LIBTORCH+x}" = x ]; then
    echo "runtime variables leaked"
    exit 9
fi
if [ "$1" = "--version" ]; then
    printf 'fibertools-rs v0.12.0\tgit-commit c881786\n'
    exit 0
fi
if [ "$1" = "predict-m6a" ] && [ "$2" = "--help" ]; then
    echo "predict help"
    exit 0
fi
exit 2
""",
            encoding="utf-8",
        )
        binary.chmod(binary.stat().st_mode | stat.S_IXUSR)

        with environment(LD_LIBRARY_PATH="/broken", LIBTORCH="/broken"):
            result = cli._validate_ft_predict_config(
                {
                    "ft_predict_backend": "external",
                    "ft_predict_executable": "bin/ft",
                    "ft_predict_expected_version": "v0.12.0",
                },
                self.root,
            )

        self.assertEqual(
            result,
            {
                "backend": "external",
                "executable": str(binary.resolve()),
                "expected_version": "0.12.0",
                "version": "0.12.0",
            },
        )

    def test_external_ft_predict_preflight_rejects_version_mismatch(self):
        binary = self.root / "ft"
        binary.write_text(
            """#!/bin/sh
if [ "$1" = "--version" ]; then
    echo "fibertools-rs v0.12.0"
    exit 0
fi
if [ "$1" = "predict-m6a" ] && [ "$2" = "--help" ]; then
    exit 0
fi
exit 2
""",
            encoding="utf-8",
        )
        binary.chmod(binary.stat().st_mode | stat.S_IXUSR)

        with self.assertRaisesRegex(
            cli.LongpipeError,
            "version mismatch: expected 0.9.0, found 0.12.0",
        ):
            cli._validate_ft_predict_config(
                {
                    "ft_predict_backend": "external",
                    "ft_predict_executable": str(binary),
                    "ft_predict_expected_version": "0.9.0",
                },
                self.root,
            )

    def test_external_ft_predict_preflight_reports_loader_failure(self):
        binary = self.root / "broken-ft"
        binary.write_text(
            """#!/bin/sh
echo "error while loading shared libraries: libtorch_cuda.so: not found"
exit 127
""",
            encoding="utf-8",
        )
        binary.chmod(binary.stat().st_mode | stat.S_IXUSR)

        with self.assertRaisesRegex(
            cli.LongpipeError,
            "libtorch_cuda.so: not found",
        ):
            cli._validate_ft_predict_config(
                {
                    "ft_predict_backend": "external",
                    "ft_predict_executable": str(binary),
                },
                self.root,
            )

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
        with (run_directory / "runsheet.csv").open(newline="") as handle:
            reader = csv.DictReader(handle)
            rows = list(reader)
        self.assertEqual(reader.fieldnames, list(cli.SETUP_RUNSHEET_COLUMNS))
        self.assertEqual(rows[0]["input_stage"], "raw")
        self.assertEqual(rows[0]["fiberseq_input"], "auto")
        self.assertIn("notes", rows[0])

        with working_directory(self.work_one):
            result, _, error = self.invoke(
                ["setup", "test-run", "--project-root", str(self.project)]
            )
        self.assertEqual(result, 2)
        self.assertIn("already exists", error)

    def test_setup_rejects_obsolete_packaged_runsheet(self):
        workflow_root = self.root / "obsolete-workflow"
        config_directory = workflow_root / "config"
        config_directory.mkdir(parents=True)
        (config_directory / "config.yaml").write_text(
            "runsheet: runsheet.csv\n", encoding="utf-8"
        )
        (config_directory / "runsheet.csv").write_text(
            "core_name,lab_name,run,cell,barcode,hifi_bam,sex,hia5\n",
            encoding="utf-8",
        )

        with (
            mock.patch.object(cli, "_workflow_root", return_value=workflow_root),
            working_directory(self.work_one),
        ):
            result, _, error = self.invoke(
                ["setup", "obsolete", "--project-root", str(self.project)]
            )

        self.assertEqual(result, 2)
        self.assertIn("packaged runsheet template is obsolete", error)
        self.assertFalse((self.work_one / "obsolete").exists())

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
        status = json.loads(
            (run_directory / "run-metadata" / "status.json").read_text()
        )
        self.assertEqual(status["command"][1], "all")
        self.assertLess(
            status["command"].index("all"),
            status["command"].index("--configfile"),
        )
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
        sample_publication = json.loads(
            (
                run_directory
                / "run-metadata"
                / "sample-publication"
                / f"{sample}.json"
            ).read_text()
        )
        self.assertEqual(sample_publication["status"], "complete")
        self.assertEqual(sample_publication["mode"], "full")
        self.assertTrue(
            (
                run_directory
                / "run-metadata"
                / "fiberseq-input"
                / f"{sample}.json"
            ).is_file()
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
            self.assertTrue(
                (invocation / "fiberseq-input" / f"{sample}.json").is_file()
            )

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

    def test_rows_runs_filtered_snapshot_and_preserves_unselected_results(self):
        run_directory = self.setup_run(self.work_one, "selected-rows")
        arguments = self.run_arguments("selected-rows")
        durable_runsheet = run_directory / "runsheet.csv"
        original_runsheet = durable_runsheet.read_bytes()

        with working_directory(self.work_one):
            result, _, error = self.invoke(arguments)
        self.assertEqual(result, 0, error)

        brain_bam = (
            run_directory / "samples" / "Brain_1" / "pbWGS" / "Brain_1.bam"
        )
        bt54_bam = run_directory / "samples" / "BT54" / "pbWGS" / "BT54.bam"
        self.assertEqual(brain_bam.read_text(), "invocation 1\n")
        self.assertEqual(bt54_bam.read_text(), "invocation 1\n")

        with working_directory(self.work_one):
            result, output, error = self.invoke(arguments + ["--rows", "2"])
        self.assertEqual(result, 0, error)
        self.assertIn("runsheet rows 2", output)
        self.assertEqual(durable_runsheet.read_bytes(), original_runsheet)
        self.assertEqual(brain_bam.read_text(), "invocation 1\n")
        self.assertEqual(bt54_bam.read_text(), "invocation 2\n")
        self.assertTrue(
            (
                run_directory
                / "run-metadata"
                / "fiberseq-input"
                / "Brain_1.json"
            ).is_file()
        )

        status = json.loads(
            (run_directory / "run-metadata" / "status.json").read_text()
        )
        self.assertEqual(status["runsheet_row_selector"], "2")
        self.assertEqual(status["runsheet_rows"], [2])
        invocation_runsheet = (
            run_directory
            / "run-metadata"
            / "invocations"
            / status["invocation_id"]
            / "runsheet.csv"
        )
        with invocation_runsheet.open(newline="") as handle:
            selected = list(csv.DictReader(handle))
        self.assertEqual(len(selected), 1)
        self.assertEqual(selected[0]["core_name"], "BT54")

    def test_rows_rejects_out_of_range_before_creating_scratch(self):
        self.setup_run(self.work_one, "bad-rows")
        with working_directory(self.work_one):
            result, _, error = self.invoke(
                self.run_arguments("bad-rows") + ["--rows", "3"]
            )
        self.assertEqual(result, 2)
        self.assertIn("has only 2 data rows", error)
        self.assertFalse((self.scratch / "bad-rows").exists())

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
        self.assertEqual(status["command"][1], "fire_only")
        self.assertLess(
            status["command"].index("fire_only"),
            status["command"].index("--configfile"),
        )
        self.assertEqual(status["status"], "success")
        fire_publication = json.loads(
            (
                run_directory
                / "run-metadata"
                / "fire-publication"
                / f"{sample}.json"
            ).read_text()
        )
        self.assertEqual(fire_publication["status"], "complete")
        self.assertEqual(fire_publication["mode"], "fire_only")

    def test_reanalysis_validation_metadata_is_published(self):
        run_directory = self.setup_run(self.work_one, "reanalysis-run")
        runsheet = run_directory / "runsheet.csv"
        runsheet.write_text(
            runsheet.read_text().replace(
                "BT53,Brain_1,r84265_20260410_231712,1_A01,bc2081,,raw,,",
                "BT53,Brain_1,r84265_20260410_231712,1_A01,bc2081,"
                "/path/to/aligned.bam,aligned,,",
            )
        )

        with working_directory(self.work_one):
            result, _, error = self.invoke(
                self.run_arguments("reanalysis-run")
            )

        self.assertEqual(result, 0, error)
        report = (
            run_directory
            / "run-metadata"
            / "reanalysis-input"
            / "Brain_1.json"
        )
        self.assertEqual(
            json.loads(report.read_text())["input_stage"], "aligned"
        )
        invocations = sorted(
            (run_directory / "run-metadata" / "invocations").iterdir()
        )
        self.assertTrue(
            (invocations[-1] / "reanalysis-input" / "Brain_1.json").is_file()
        )

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
        self.assertTrue(
            (
                run_directory
                / "run-metadata"
                / "fiberseq-input"
                / "Brain_1.json"
            ).is_file()
        )
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
        self.assertFalse(
            (
                run_directory
                / "run-metadata"
                / "fiberseq-input"
                / "Brain_1.json"
            ).exists()
        )
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

    def test_completed_sample_publishes_when_another_sample_fails(self):
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
        self.assertEqual(result, 1)
        self.assertIn("Run failed with exit code 1", error)
        self.assertEqual(pileup.read_text(), "invocation 2\n")
        self.assertEqual(
            (
                run_directory
                / "samples"
                / "BT54"
                / "pbWGS"
                / "BT54.bam"
            ).read_text(),
            "invocation 1\n",
        )
        status = json.loads(
            (run_directory / "run-metadata" / "status.json").read_text()
        )
        self.assertEqual(status["status"], "failed")

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
