from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import yaml


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


class ConfigSchemaTest(unittest.TestCase):
    def test_fire_minimum_mean_depth_defaults_to_eight(self):
        schema = yaml.safe_load(
            (REPOSITORY_ROOT / "schemas" / "config.schema.yaml").read_text()
        )

        self.assertEqual(
            schema["properties"]["fire_min_mean_depth"]["default"], 8.0
        )

    def test_fire_version_defaults_to_shared_release(self):
        schema = yaml.safe_load(
            (REPOSITORY_ROOT / "schemas" / "config.schema.yaml").read_text()
        )

        self.assertEqual(
            schema["properties"]["fire_version"]["default"], "v0.1.2"
        )

    def test_shared_miniwdl_deployment_is_pinned(self):
        environment = yaml.safe_load(
            (REPOSITORY_ROOT / "deploy" / "miniwdl" / "environment.yaml")
            .read_text(encoding="utf-8")
        )
        dependencies = environment["dependencies"]
        pip_dependencies = next(
            dependency["pip"]
            for dependency in dependencies
            if isinstance(dependency, dict) and "pip" in dependency
        )
        config = yaml.safe_load(
            (REPOSITORY_ROOT / "config" / "config.yaml").read_text(
                encoding="utf-8"
            )
        )

        self.assertIn("miniwdl=1.15.0", dependencies)
        self.assertIn("miniwdl-slurm==0.4.0", pip_dependencies)
        self.assertEqual(
            config["miniwdl_executable"],
            "/fh/fast/furlan_s/grp/software/miniwdl/bin/miniwdl",
        )

    def test_pacbio_references_use_v3_bundle(self):
        config = yaml.safe_load(
            (REPOSITORY_ROOT / "config" / "config.yaml").read_text(
                encoding="utf-8"
            )
        )
        reference_keys = (
            "ref_fasta",
            "fire_ref",
            "chrom_sizes",
            "ref_map_file",
            "tertiary_map_file",
            "cnv_excluded_regions",
            "expected_cn_xx",
            "expected_cn_xy",
        )

        for key in reference_keys:
            with self.subTest(key=key):
                self.assertIn("/hifi-wdl-resources-v3.1.0/", config[key])
        self.assertTrue(
            config["ref_map_file"].endswith("GRCh38.ref_map.v3p1p0.tsv")
        )
        self.assertTrue(
            config["tertiary_map_file"].endswith(
                "GRCh38.tertiary_map.v3p1p0.tsv"
            )
        )

    def test_slurm_profile_allows_ten_wdl_drivers(self):
        profile = yaml.safe_load(
            (REPOSITORY_ROOT / "profiles" / "slurm" / "config.yaml").read_text()
        )

        self.assertEqual(profile["resources"]["wdl_driver"], 10)
        self.assertGreater(
            profile["cores"], profile["resources"]["wdl_driver"]
        )
        self.assertEqual(profile["local-cores"], profile["cores"])

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

    def test_fire_checkpoint_preserves_and_publishes_pileup(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (
                run_directory,
                scratch_directory,
                inputs,
                config_path,
                environment,
            ) = self._fixture(
                root,
                "fire-run",
                "yes",
                pixi_executable="inputs/fake-pixi",
            )

            Path(environment["LONGPIPE_RUNSHEET"]).write_text(
                "core_name,lab_name,run,cell,barcode,hifi_bam,sex,hia5\n"
                "sample,,,,,inputs/sample.bam,F,yes\n"
                "current_key,,,,,inputs/sample.bam,F,yes\n",
                encoding="utf-8",
            )

            fake_pixi = inputs / "fake-pixi"
            fake_pixi.write_text(
                "#!/usr/bin/env bash\n"
                "set -euo pipefail\n"
                "(\n"
                "IFS= read -r _\n"
                "while IFS=$'\\t' read -r sample _; do\n"
                "  mkdir -p \"results/$sample\"\n"
                "  printf 'pileup\\n' > \"results/$sample/"
                "$sample-fire-v0.1.1-pileup.bed.gz\"\n"
                "  printf 'index\\n' > \"results/$sample/"
                "$sample-fire-v0.1.1-pileup.bed.gz.tbi\"\n"
                "done\n"
                ") < config.tbl\n",
                encoding="utf-8",
            )
            fake_pixi.chmod(0o755)

            depth_outputs = {
                "sample": {"humanwgs_singleton.stat_mean_depth": 20.0},
                "current_key": {
                    "humanwgs_singleton.stat_depth_mean": 21.0
                },
            }
            for sample, outputs in depth_outputs.items():
                pbwgs = scratch_directory / "pbWGS" / sample
                pbwgs.mkdir(parents=True)
                (pbwgs / "outputs.json").write_text(
                    json.dumps({"outputs": outputs}),
                    encoding="utf-8",
                )
                wgs = scratch_directory / "results" / "wgs" / sample
                wgs.mkdir(parents=True)
                (wgs / f"{sample}.GRCh38.haplotagged.bam").touch()

            target = (
                "results/fire/pileups/"
                "sample-fire-v0.1.1-pileup.bed.gz"
            )
            completed = subprocess.run(
                [
                    "snakemake",
                    target,
                    "--snakefile",
                    str(REPOSITORY_ROOT / "Snakefile"),
                    "--configfile",
                    str(config_path),
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

            fire_config = yaml.safe_load(
                (scratch_directory / "fibertools" / "FIRE" / "config.yaml")
                .read_text(encoding="utf-8")
            )
            self.assertIs(fire_config["full-version"], True)

            resumed = subprocess.run(
                [
                    "snakemake",
                    target,
                    "--snakefile",
                    str(REPOSITORY_ROOT / "Snakefile"),
                    "--configfile",
                    str(config_path),
                    "--cores",
                    "1",
                    "--forcerun",
                    "fire",
                ],
                cwd=run_directory,
                env=environment,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            self.assertEqual(resumed.returncode, 0, resumed.stdout)

            raw_pileup = (
                scratch_directory
                / "fibertools"
                / "FIRE"
                / "results"
                / "sample"
                / "sample-fire-v0.1.1-pileup.bed.gz"
            )
            self.assertEqual(raw_pileup.read_text(), "pileup\n")
            self.assertEqual(
                (scratch_directory / target).read_text(), "pileup\n"
            )
            staged_pileup = (
                scratch_directory
                / "intermediate"
                / "fire"
                / "pileups"
                / "sample-fire-v0.1.1-pileup.bed.gz"
            )
            self.assertEqual(staged_pileup.read_text(), "pileup\n")
            self.assertEqual(
                raw_pileup.stat().st_ino,
                staged_pileup.stat().st_ino,
            )
            raw_index = raw_pileup.with_name(raw_pileup.name + ".tbi")
            staged_index = staged_pileup.with_name(staged_pileup.name + ".tbi")
            self.assertEqual(raw_index.read_text(), "index\n")
            self.assertEqual(raw_index.stat().st_ino, staged_index.stat().st_ino)

    def _fixture(
        self,
        root: Path,
        run_name: str,
        hia5: str,
        pixi_executable: str = "pixi",
    ) -> tuple[Path, Path, Path, Path, dict[str, str]]:
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
            "pixi_executable": pixi_executable,
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
        return (
            run_directory,
            scratch_directory,
            inputs,
            config_path,
            environment,
        )

    def _dry_run(self, root: Path, run_name: str, hia5: str) -> None:
        (
            run_directory,
            scratch_directory,
            _,
            config_path,
            environment,
        ) = self._fixture(root, run_name, hia5)
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
        self.assertIn("checkpoint fire_eligibility", completed.stdout)
        if hia5 == "yes":
            self.assertIn("rule add_nucleosomes", completed.stdout)
        else:
            self.assertNotIn("rule add_nucleosomes", completed.stdout)
