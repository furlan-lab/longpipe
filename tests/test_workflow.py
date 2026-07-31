from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import yaml


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


class ConfigSchemaTest(unittest.TestCase):
    def test_ccre_peaks_requires_nonempty_safe_name_mapping(self):
        from jsonschema import Draft202012Validator

        schema = yaml.safe_load(
            (REPOSITORY_ROOT / "schemas" / "config.schema.yaml").read_text()
        )
        validator = Draft202012Validator(
            schema["properties"]["ccre_peaks"]
        )

        self.assertFalse(
            list(
                validator.iter_errors(
                    {
                        "encode": "refs/encode.bed.gz",
                        "brain_v2": "/refs/brain.bed.gz",
                    }
                )
            )
        )
        self.assertTrue(list(validator.iter_errors("refs/encode.bed.gz")))
        self.assertTrue(list(validator.iter_errors({})))
        self.assertTrue(
            list(validator.iter_errors({"invalid/name": "refs/encode.bed.gz"}))
        )


@unittest.skipUnless(shutil.which("snakemake"), "snakemake is not installed")
class WorkflowDagTest(unittest.TestCase):
    def test_read_and_fiber_branches_build_in_scratch(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for run_name, hia5 in (("read-run", "no"), ("fiber-run", "yes")):
                with self.subTest(hia5=hia5):
                    self._dry_run(root, run_name, hia5)

    def _dry_run(self, root: Path, run_name: str, hia5: str) -> None:
        run_directory = root / "project" / run_name
        scratch_directory = root / "scratch" / run_name
        inputs = run_directory / "inputs"
        invocation = scratch_directory / ".longpipe" / "invocations" / "test"
        inputs.mkdir(parents=True)
        invocation.mkdir(parents=True)

        input_names = [
            "sample.bam",
            "ref.fa",
            "ref.fa.fai",
            "hla.fa",
            "chrom.sizes",
            "encode-ccre.bed.gz",
            "brain-ccre.bed.gz",
            "ref-map.tsv",
            "tertiary-map.tsv",
            "workflow.wdl",
            "fire.toml",
            "cnv.bed.gz",
            "xx.bed",
            "xy.bed",
            "copy_pbWDL",
            "shrink_fire_bed",
        ]
        for name in input_names:
            (inputs / name).touch()

        runsheet = invocation / "runsheet.csv"
        runsheet.write_text(
            "core_name,lab_name,run,cell,barcode,hifi_bam,sex,hia5\n"
            f"sample,,,,,inputs/sample.bam,F,{hia5}\n",
            encoding="utf-8",
        )
        config = {
            "rawdir": "inputs",
            "runsheet": "runsheet.csv",
            "default_sex": "FEMALE",
            "hia5_gates_fiberseq": True,
            "hla_mode": "reads",
            "scatac_sorted": True,
            "fire_version": "v0.1.1",
            "wdl_preemptible": True,
            "ref_fasta": "inputs/ref.fa",
            "fire_ref": "inputs/ref.fa",
            "hla_ref_fasta": "inputs/hla.fa",
            "chrom_sizes": "inputs/chrom.sizes",
            "ccre_peaks": {
                "encode": "inputs/encode-ccre.bed.gz",
                "brain": "inputs/brain-ccre.bed.gz",
            },
            "ref_map_file": "inputs/ref-map.tsv",
            "tertiary_map_file": "inputs/tertiary-map.tsv",
            "wdl_singleton": "inputs/workflow.wdl",
            "fire_pixi": "inputs/fire.toml",
            "cnv_excluded_regions": "inputs/cnv.bed.gz",
            "expected_cn_xx": "inputs/xx.bed",
            "expected_cn_xy": "inputs/xy.bed",
            "miniwdl_executable": "miniwdl",
            "pixi_executable": "pixi",
            "copy_pbwdl_tool": "inputs/copy_pbWDL",
            "shrink_fire_bed_tool": "inputs/shrink_fire_bed",
        }
        config_path = invocation / "config.yaml"
        config_path.write_text(
            yaml.safe_dump(config, sort_keys=False),
            encoding="utf-8",
        )

        environment = os.environ.copy()
        environment.update(
            {
                "LONGPIPE_RUN_DIR": str(run_directory),
                "LONGPIPE_RUN_NAME": run_name,
                "LONGPIPE_RUNSHEET": str(runsheet),
                "LONGPIPE_SCRATCH_DIR": str(scratch_directory),
            }
        )
        completed = subprocess.run(
            [
                "snakemake",
                "--snakefile",
                str(REPOSITORY_ROOT / "Snakefile"),
                "--configfile",
                str(config_path),
                "--dry-run",
                "--cores",
                "1",
            ],
            cwd=run_directory,
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        self.assertEqual(completed.returncode, 0, completed.stdout)
        self.assertTrue((scratch_directory / ".snakemake").is_dir())
        self.assertFalse((run_directory / ".snakemake").exists())
        if hia5 == "yes":
            for peak_set in ("encode", "brain"):
                expected = (
                    f"results/fire/ccre/{peak_set}/"
                    "sample_scATAC_max_scores.bed"
                )
                self.assertIn(expected, completed.stdout)
