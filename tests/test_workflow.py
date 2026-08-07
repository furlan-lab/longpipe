from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
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

    def test_haplotype_analysis_defaults_are_explicit(self):
        schema = yaml.safe_load(
            (REPOSITORY_ROOT / "schemas" / "config.schema.yaml").read_text()
        )
        properties = schema["properties"]

        self.assertEqual(
            properties["fire_haplotype_min_coverage"]["default"], 10
        )
        self.assertEqual(
            properties["fire_haplotype_max_fdr"]["default"], 0.05
        )
        self.assertEqual(properties["fire_haplotype_qvalue"]["default"], 0.05)

    def test_native_fire_source_is_pinned_and_licensed(self):
        schema = yaml.safe_load(
            (REPOSITORY_ROOT / "schemas" / "config.schema.yaml").read_text()
        )
        expected_commit = "e7f2aae1f59bbb74c46e1f27936da380cb074704"

        self.assertEqual(
            schema["properties"]["fire_upstream_commit"]["const"],
            expected_commit,
        )
        license_text = (
            REPOSITORY_ROOT / "workflow" / "fire" / "LICENSE"
        ).read_text(encoding="utf-8")
        self.assertIn("MIT License", license_text)
        snakefile = (REPOSITORY_ROOT / "Snakefile").read_text(encoding="utf-8")
        self.assertNotIn("fire_pixi", snakefile)
        self.assertNotIn("pixi_executable", schema["properties"])

    def test_fire_excludes_require_a_nonempty_path_list(self):
        from jsonschema import Draft202012Validator

        schema = yaml.safe_load(
            (REPOSITORY_ROOT / "schemas" / "config.schema.yaml").read_text()
        )
        validator = Draft202012Validator(
            schema["properties"]["fire_exclude_regions"]
        )

        self.assertFalse(
            list(validator.iter_errors(["gap.bed.gz", "blacklist.bed.gz"]))
        )
        self.assertTrue(list(validator.iter_errors([])))
        self.assertTrue(list(validator.iter_errors("gap.bed.gz")))

    def test_fire_python_environment_preserves_numeric_pins(self):
        environment = yaml.safe_load(
            (REPOSITORY_ROOT / "envs" / "fire-python.yaml").read_text(
                encoding="utf-8"
            )
        )
        dependencies = environment["dependencies"]
        pip_dependencies = next(
            dependency["pip"]
            for dependency in dependencies
            if isinstance(dependency, dict) and "pip" in dependency
        )

        self.assertIn("numpy=1.24", dependencies)
        self.assertIn("pandas=1.4", dependencies)
        self.assertEqual(pip_dependencies, ["polars-lts-cpu==1.7.1"])

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
        self.assertNotIn("fire_driver", profile["resources"])
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


class FirePileupComparisonTest(unittest.TestCase):
    def test_comparison_uses_decompressed_content(self):
        import gzip

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            legacy = root / "legacy.bed.gz"
            native = root / "native.bed.gz"
            content = b"#chrom\tstart\tend\nchr1\t0\t10\n"
            with gzip.open(legacy, "wb", compresslevel=1) as handle:
                handle.write(content)
            with gzip.open(native, "wb", compresslevel=9) as handle:
                handle.write(content)

            completed = subprocess.run(
                [
                    sys.executable,
                    str(
                        REPOSITORY_ROOT
                        / "workflow"
                        / "scripts"
                        / "compare_fire_pileups.py"
                    ),
                    str(legacy),
                    str(native),
                ],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )

            self.assertEqual(completed.returncode, 0, completed.stdout)
            self.assertIn("identical after decompression", completed.stdout)


@unittest.skipUnless(shutil.which("snakemake"), "snakemake is not installed")
class WorkflowDagTest(unittest.TestCase):
    def test_read_and_fiber_branches_build_in_scratch(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for run_name, hia5 in (("read-run", "no"), ("fiber-run", "yes")):
                with self.subTest(hia5=hia5):
                    self._dry_run(root, run_name, hia5)

    def test_fire_only_uses_native_dag(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (
                run_directory,
                scratch_directory,
                _,
                config_path,
                environment,
            ) = self._fixture(root, "fire-run", "yes")

            pbwgs = scratch_directory / "pbWGS" / "sample"
            pbwgs.mkdir(parents=True)
            (pbwgs / "outputs.json").write_text(
                json.dumps(
                    {
                        "outputs": {
                            "humanwgs_singleton.stat_depth_mean": 20.0
                        }
                    }
                ),
                encoding="utf-8",
            )
            wgs = scratch_directory / "results" / "wgs" / "sample"
            wgs.mkdir(parents=True)
            (wgs / "sample.GRCh38.haplotagged.bam").touch()
            (wgs / "sample.GRCh38.haplotagged.bam.bai").touch()

            eligibility_directory = scratch_directory / "intermediate" / "fire"
            eligibility_directory.mkdir(parents=True)
            (eligibility_directory / "eligibility.json").write_text(
                json.dumps(
                    {
                        "minimum_mean_depth": 8.0,
                        "eligible_samples": ["sample"],
                        "samples": [],
                    }
                ),
                encoding="utf-8",
            )
            report = scratch_directory / "results" / "fire" / "eligibility.tsv"
            report.parent.mkdir(parents=True)
            report.write_text(
                "sample\tmean_depth\tminimum_mean_depth\teligible\treason\n"
                "sample\t20\t8\tyes\teligible\n",
                encoding="utf-8",
            )

            completed = subprocess.run(
                [
                    "snakemake",
                    "fire_only",
                    "--snakefile",
                    str(REPOSITORY_ROOT / "Snakefile"),
                    "--configfile",
                    str(config_path),
                    "--dry-run",
                    "--cores",
                    "11",
                ],
                cwd=run_directory,
                env=environment,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            self.assertEqual(completed.returncode, 0, completed.stdout)
            native_rules = (
                "fire_filtered_cram",
                "fire_genome_bedgraph",
                "fire_coverage",
                "fire_fiber_locations_chromosome",
                "fire_fiber_locations",
                "fire_exclude_from_shuffle",
                "fire_shuffled_fiber_locations_chromosome",
                "fire_shuffled_pileup_chromosome",
                "fire_shuffled_pileup",
                "fire_fdr_table",
                "fire_pileup_chromosome",
                "fire_fdr_track_chromosome",
                "fire_pileup",
                "fire_bigwig",
                "fire_haplotype_tracks",
                "fire_haplotype_differences",
            )
            for rule in native_rules:
                with self.subTest(rule=rule):
                    self.assertIn(f"rule {rule}", completed.stdout)
            self.assertNotIn("pixi", completed.stdout.lower())
            self.assertNotIn("rule scatac_map", completed.stdout)
            self.assertNotIn("rule sawfish", completed.stdout)
            self.assertNotIn("rule hifihla", completed.stdout)

    def _fixture(
        self,
        root: Path,
        run_name: str,
        hia5: str,
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
            "gap.bed.gz",
            "blacklist.bed.gz",
            "segmental-duplications.bed.gz",
            "cnv.bed.gz",
            "xx.bed",
            "xy.bed",
            "copy_pbWDL",
            "shrink_fire_bed",
        ]
        for name in input_names:
            (inputs / name).touch()
        (inputs / "ref.fa.fai").write_text(
            "chr1\t100000\t0\t80\t81\n",
            encoding="utf-8",
        )
        (inputs / "chrom.sizes").write_text(
            "chr1\t100000\n",
            encoding="utf-8",
        )

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
            "fire_version": "v0.1.2",
            "fire_upstream_commit": (
                "e7f2aae1f59bbb74c46e1f27936da380cb074704"
            ),
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
            "fire_exclude_regions": [
                "inputs/gap.bed.gz",
                "inputs/blacklist.bed.gz",
                "inputs/segmental-duplications.bed.gz",
            ],
            "cnv_excluded_regions": "inputs/cnv.bed.gz",
            "expected_cn_xx": "inputs/xx.bed",
            "expected_cn_xy": "inputs/xy.bed",
            "miniwdl_executable": "miniwdl",
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
