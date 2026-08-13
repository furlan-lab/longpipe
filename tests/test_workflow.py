from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import tomllib
import unittest
from pathlib import Path

import yaml


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


class ConfigSchemaTest(unittest.TestCase):
    FIBERTOOLS_VERSION = "0.12.0"

    def test_release_version_is_0_2_5(self):
        project = tomllib.loads(
            (REPOSITORY_ROOT / "pyproject.toml").read_text(encoding="utf-8")
        )

        self.assertEqual(project["project"]["version"], "0.2.5")

    def test_cached_helper_scripts_are_rule_inputs(self):
        workflow_sources = "\n".join(
            path.read_text(encoding="utf-8")
            for path in (
                REPOSITORY_ROOT / "Snakefile",
                REPOSITORY_ROOT / "workflow" / "rules" / "fire.smk",
            )
        )

        self.assertNotIn("{params.script", workflow_sources)
        self.assertEqual(workflow_sources.count("{input.script:q}"), 7)
        self.assertEqual(workflow_sources.count("script_sha256="), 7)
        self.assertIn("return ancient(script)", workflow_sources)

    def test_runsheet_accepts_free_text_notes(self):
        from jsonschema import Draft202012Validator

        schema = yaml.safe_load(
            (REPOSITORY_ROOT / "schemas" / "runsheet.schema.yaml").read_text()
        )
        validator = Draft202012Validator(schema["properties"]["notes"])

        self.assertFalse(list(validator.iter_errors("healthy donor; batch 2")))
        self.assertFalse(list(validator.iter_errors("")))
        self.assertTrue(list(validator.iter_errors(42)))

    def test_ft_predict_defaults_to_path_lookup(self):
        schema = yaml.safe_load(
            (REPOSITORY_ROOT / "schemas" / "config.schema.yaml").read_text()
        )

        self.assertEqual(
            schema["properties"]["ft_predict_executable"]["default"], "ft"
        )
        self.assertEqual(
            schema["properties"]["ft_predict_backend"]["enum"],
            ["conda", "external"],
        )

    def test_runsheet_accepts_explicit_fiberseq_input_types(self):
        from jsonschema import Draft202012Validator

        schema = yaml.safe_load(
            (REPOSITORY_ROOT / "schemas" / "runsheet.schema.yaml").read_text()
        )
        validator = Draft202012Validator(
            schema["properties"]["fiberseq_input"]
        )

        self.assertFalse(list(validator.iter_errors("m6a")))
        self.assertFalse(list(validator.iter_errors("kinetics")))
        self.assertFalse(list(validator.iter_errors("auto")))
        self.assertFalse(list(validator.iter_errors("")))
        self.assertTrue(list(validator.iter_errors("unknown")))

    def test_runsheet_accepts_explicit_input_stages(self):
        from jsonschema import Draft202012Validator

        schema = yaml.safe_load(
            (REPOSITORY_ROOT / "schemas" / "runsheet.schema.yaml").read_text()
        )
        validator = Draft202012Validator(schema["properties"]["input_stage"])

        for value in ("", "raw", "aligned", "pbwgs", "PBWGS"):
            with self.subTest(value=value):
                self.assertFalse(list(validator.iter_errors(value)))
        self.assertTrue(list(validator.iter_errors("unknown")))

    def test_fiberseq_detection_defaults_match_ft_validate(self):
        schema = yaml.safe_load(
            (REPOSITORY_ROOT / "schemas" / "config.schema.yaml").read_text()
        )
        properties = schema["properties"]

        self.assertEqual(
            properties["fiberseq_detection_reads"]["default"],
            5000,
        )
        self.assertEqual(
            properties["fiberseq_detection_min_fraction"]["default"],
            0.5,
        )

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

    def test_fibertools_conda_environments_share_version(self):
        fibertools_dependencies = {}
        for environment_path in sorted((REPOSITORY_ROOT / "envs").glob("*.yaml")):
            environment = yaml.safe_load(environment_path.read_text(encoding="utf-8"))
            pins = [
                dependency
                for dependency in environment["dependencies"]
                if isinstance(dependency, str)
                and dependency.startswith("fibertools-rs=")
            ]
            if pins:
                fibertools_dependencies[environment_path.name] = pins

        self.assertEqual(
            fibertools_dependencies,
            {
                "fibertools.yaml": [
                    f"fibertools-rs={self.FIBERTOOLS_VERSION}"
                ],
                "fire-core.yaml": [
                    f"fibertools-rs={self.FIBERTOOLS_VERSION}"
                ],
            },
        )

        fire_core = yaml.safe_load(
            (REPOSITORY_ROOT / "envs" / "fire-core.yaml").read_text(
                encoding="utf-8"
            )
        )["dependencies"]
        self.assertIn("samtools=1.23.1", fire_core)
        self.assertIn("htslib=1.23.1", fire_core)
        self.assertIn("mosdepth=0.3.7", fire_core)

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


class FiberseqInputDetectionTest(unittest.TestCase):
    def test_auto_selects_each_supported_action(self):
        cases = {
            "m6a": (
                ["MM:Z:A+a,0;", "ML:B:C,200"],
                "add_nucleosomes",
            ),
            "m6a_with_retained_kinetics": (
                [
                    "MM:Z:A+a,0;",
                    "ML:B:C,200",
                    "fi:B:C,1,1,1,1",
                    "ri:B:C,1,1,1,1",
                    "fp:B:C,1,1,1,1",
                    "rp:B:C,1,1,1,1",
                ],
                "add_nucleosomes",
            ),
            "complete": (
                [
                    "MM:Z:A+a,0;",
                    "ML:B:C,200",
                    "ns:B:I,0",
                    "nl:B:I,100",
                    "as:B:I,100",
                    "al:B:I,50",
                ],
                "use_existing",
            ),
            "complete_ma": (
                [
                    "MM:Z:A+a,0;",
                    "ML:B:C,200",
                    "MA:Z:4;nuc+:1-2;msp+:3-2",
                ],
                "use_existing",
            ),
            "unrelated_ma": (
                [
                    "MM:Z:A+a,0;",
                    "ML:B:C,200",
                    "MA:Z:4;ctcf+:1-2",
                ],
                "add_nucleosomes",
            ),
            "kinetics": (
                [
                    "fi:B:C,1,1,1,1",
                    "ri:B:C,1,1,1,1",
                    "fp:B:C,1,1,1,1",
                    "rp:B:C,1,1,1,1",
                ],
                "predict_m6a",
            ),
            "cpg_with_kinetics": (
                [
                    "MM:Z:C+m,0;",
                    "ML:B:C,200",
                    "fi:B:C,1,1,1,1",
                    "ri:B:C,1,1,1,1",
                    "fp:B:C,1,1,1,1",
                    "rp:B:C,1,1,1,1",
                ],
                "predict_m6a",
            ),
        }
        for name, (tags, expected) in cases.items():
            with self.subTest(name=name):
                completed, report = self._detect(tags)
                self.assertEqual(completed.returncode, 0, completed.stdout)
                self.assertEqual(report["selected_action"], expected)
                self.assertEqual(report["reads_tested"], 1)

    def test_ma_detection_requires_both_nucleosome_and_msp_types(self):
        completed, report = self._detect(
            [
                "MM:Z:A+a,0;",
                "ML:B:C,200",
                "MA:Z:4;nuc+:1-2",
            ]
        )

        self.assertEqual(completed.returncode, 0, completed.stdout)
        self.assertEqual(report["selected_action"], "add_nucleosomes")
        self.assertEqual(report["ma_annotation_reads"], 0)

    def test_auto_rejects_bam_without_m6a_or_kinetics(self):
        completed, report = self._detect(["rq:f:0.999"])

        self.assertEqual(completed.returncode, 2, completed.stdout)
        self.assertIsNone(report)
        self.assertIn("neither sufficient A+a MM/ML", completed.stdout)

    def test_explicit_override_still_validates_required_tags(self):
        completed, report = self._detect(
            ["MM:Z:A+a,0;", "ML:B:C,200"],
            mode="kinetics",
        )

        self.assertEqual(completed.returncode, 2, completed.stdout)
        self.assertIsNone(report)
        self.assertIn("requires fi, ri, fp, and rp", completed.stdout)

    def _detect(
        self,
        tags: list[str],
        mode: str = "auto",
    ) -> tuple[subprocess.CompletedProcess[str], dict | None]:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            alignments = root / "sample.sam"
            alignment = "\t".join(
                [
                    "read1",
                    "4",
                    "*",
                    "0",
                    "0",
                    "*",
                    "*",
                    "0",
                    "0",
                    "AAAA",
                    "IIII",
                    *tags,
                ]
            )
            alignments.write_text(f"{alignment}\n", encoding="utf-8")
            fake_samtools = root / "samtools"
            fake_samtools.write_text(
                "#!/usr/bin/env python3\n"
                "import sys\n"
                "from pathlib import Path\n"
                "sys.stdout.write(Path(sys.argv[-1]).read_text())\n",
                encoding="utf-8",
            )
            fake_samtools.chmod(0o755)
            output = root / "detection.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(
                        REPOSITORY_ROOT
                        / "workflow"
                        / "scripts"
                        / "detect_fiberseq_input.py"
                    ),
                    "--input",
                    str(alignments),
                    "--output",
                    str(output),
                    "--sample",
                    "sample",
                    "--mode",
                    mode,
                    "--samtools",
                    str(fake_samtools),
                ],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            report = json.loads(output.read_text()) if output.is_file() else None
            return completed, report


class ReanalysisInputValidationTest(unittest.TestCase):
    def test_valid_pbwgs_bam_records_version_and_action(self):
        completed, report = self._validate("pbwgs")

        self.assertEqual(completed.returncode, 0, completed.stdout)
        self.assertEqual(report["source_workflow_version"], "v2.0.7")
        self.assertEqual(report["action"], "reuse_pbwgs_outputs")
        self.assertTrue(report["reference_dictionary_matches"])
        self.assertEqual(report["sampled_records"]["mapped_primary"], 1)

    def test_aligned_bam_is_marked_for_current_wdl_realign(self):
        completed, report = self._validate("aligned")

        self.assertEqual(completed.returncode, 0, completed.stdout)
        self.assertIsNone(report["source_workflow_version"])
        self.assertEqual(report["action"], "realign_with_current_pbwgs")

    def test_reference_dictionary_mismatch_fails(self):
        completed, report = self._validate("aligned", reference_length=99)

        self.assertEqual(completed.returncode, 2, completed.stdout)
        self.assertIsNone(report)
        self.assertIn("reference dictionary does not match", completed.stdout)

    def test_reference_dictionary_order_mismatch_fails(self):
        completed, report = self._validate(
            "aligned",
            reference_contigs=(("chr2", 200), ("chr1", 100)),
        )

        self.assertEqual(completed.returncode, 2, completed.stdout)
        self.assertIsNone(report)
        self.assertIn("order_mismatch=True", completed.stdout)

    def _validate(
        self,
        stage: str,
        reference_length: int = 100,
        reference_contigs: tuple[tuple[str, int], ...] | None = None,
    ) -> tuple[subprocess.CompletedProcess[str], dict | None]:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bam = root / "sample.GRCh38.haplotagged.bam"
            bam.touch()
            reference_index = root / "ref.fa.fai"
            if reference_contigs is None:
                reference_contigs = (
                    ("chr1", reference_length),
                    ("chr2", 200),
                )
            reference_index.write_text(
                "".join(
                    f"{name}\t{length}\t0\t80\t81\n"
                    for name, length in reference_contigs
                ),
                encoding="utf-8",
            )
            outputs = root / "outputs.json"
            outputs.write_text(
                json.dumps(
                    {
                        "humanwgs_singleton.workflow_version": "v2.0.7",
                        "humanwgs_singleton.merged_haplotagged_bam": str(bam),
                    }
                ),
                encoding="utf-8",
            )
            fake_samtools = root / "samtools"
            fake_samtools.write_text(
                "#!/usr/bin/env python3\n"
                "import sys\n"
                "if sys.argv[1] == 'quickcheck':\n"
                "    raise SystemExit(0)\n"
                "if '-H' in sys.argv:\n"
                "    print('@HD\\tVN:1.6\\tSO:coordinate')\n"
                "    print('@SQ\\tSN:chr1\\tLN:100')\n"
                "    print('@SQ\\tSN:chr2\\tLN:200')\n"
                "else:\n"
                "    print('read1\\t0\\tchr1\\t1\\t60\\t4=\\t*\\t0\\t0\\tAAAA\\tIIII')\n",
                encoding="utf-8",
            )
            fake_samtools.chmod(0o755)
            report_path = root / "validation.json"
            command = [
                sys.executable,
                str(
                    REPOSITORY_ROOT
                    / "workflow"
                    / "scripts"
                    / "validate_reanalysis_input.py"
                ),
                "--input",
                str(bam),
                "--output",
                str(report_path),
                "--reference-index",
                str(reference_index),
                "--sample",
                "sample",
                "--stage",
                stage,
                "--samtools",
                str(fake_samtools),
            ]
            if stage == "pbwgs":
                command.extend(("--outputs-json", str(outputs)))
            completed = subprocess.run(
                command,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            report = (
                json.loads(report_path.read_text())
                if report_path.is_file()
                else None
            )
            return completed, report


@unittest.skipUnless(shutil.which("snakemake"), "snakemake is not installed")
class WorkflowDagTest(unittest.TestCase):
    def test_read_and_fiber_branches_build_in_scratch(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for run_name, hia5 in (("read-run", "no"), ("fiber-run", "yes")):
                with self.subTest(hia5=hia5):
                    self._dry_run(root, run_name, hia5)

    def test_cached_detector_reinstall_does_not_break_resume(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            workflow_root = root / "workflow-source"
            shutil.copytree(
                REPOSITORY_ROOT,
                workflow_root,
                ignore=shutil.ignore_patterns(".git", ".pytest_cache"),
            )
            (
                run_directory,
                scratch_directory,
                inputs,
                config_path,
                environment,
            ) = self._fixture(root, "detector-resume", "yes")
            fake_samtools = inputs / "samtools"
            fake_samtools.write_text(
                "#!/bin/sh\n"
                "printf 'read1\\t4\\t*\\t0\\t0\\t*\\t*\\t0\\t0\\tAAAA"
                "\\tIIII\\tMM:Z:A+a,0;\\tML:B:C,255\\n'\n",
                encoding="utf-8",
            )
            fake_samtools.chmod(0o755)
            environment["PATH"] = f"{inputs}{os.pathsep}{environment['PATH']}"
            target = (
                scratch_directory
                / "intermediate"
                / "fiberseq"
                / "input-detection"
                / "sample.json"
            )
            command = [
                "snakemake",
                str(target.relative_to(scratch_directory)),
                "--snakefile",
                str(workflow_root / "Snakefile"),
                "--configfile",
                str(config_path),
                "--cores",
                "1",
            ]

            first = subprocess.run(
                command,
                cwd=run_directory,
                env=environment,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            self.assertEqual(first.returncode, 0, first.stdout)
            detector = (
                workflow_root
                / "workflow"
                / "scripts"
                / "detect_fiberseq_input.py"
            )
            detector_contents = detector.read_bytes()
            detector.write_bytes(detector_contents)
            resumed = subprocess.run(
                [*command, "--dry-run"],
                cwd=run_directory,
                env=environment,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            detector.write_bytes(detector_contents + b"\n")
            changed = subprocess.run(
                [*command, "--dry-run"],
                cwd=run_directory,
                env=environment,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )

        self.assertEqual(resumed.returncode, 0, resumed.stdout)
        self.assertNotIn("checkpoint detect_fiberseq_input", resumed.stdout)
        self.assertIn("Nothing to be done", resumed.stdout)
        self.assertEqual(changed.returncode, 0, changed.stdout)
        self.assertIn("checkpoint detect_fiberseq_input", changed.stdout)
        self.assertIn("Params have changed", changed.stdout)

    def test_targeting_one_sample_publication_does_not_require_other_sample(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (
                run_directory,
                _,
                inputs,
                config_path,
                environment,
            ) = self._fixture(root, "independent-samples", "no")
            (inputs / "second.bam").touch()
            Path(environment["LONGPIPE_RUNSHEET"]).write_text(
                "core_name,lab_name,run,cell,barcode,hifi_bam,input_stage,"
                "pbwgs_outputs_json,sex,hia5,fiberseq_input\n"
                "sample,,,,,inputs/sample.bam,raw,,F,no,auto\n"
                "second,,,,,inputs/second.bam,raw,,F,no,auto\n",
                encoding="utf-8",
            )
            target = str(
                (
                    run_directory
                    / "run-metadata"
                    / "sample-publication"
                    / "sample.json"
                ).resolve()
            )
            completed = self._run_dry_run(
                run_directory,
                config_path,
                environment,
                target=target,
            )

        self.assertIn("rule publish_sample_results", completed.stdout)
        self.assertIn("wildcards: sample=sample", completed.stdout)
        self.assertNotIn("second.bam", completed.stdout)
        self.assertNotIn("wildcards: sample=second", completed.stdout)

    def test_sample_publication_rule_executes_atomically(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (
                run_directory,
                scratch_directory,
                _,
                config_path,
                environment,
            ) = self._fixture(root, "publish-one", "no")

            wgs = scratch_directory / "results" / "wgs" / "sample"
            wgs.mkdir(parents=True)
            bam = wgs / "sample.GRCh38.haplotagged.bam"
            bam.write_text("bam\n", encoding="utf-8")
            (wgs / "outputs.json").write_text(
                json.dumps({"haplotagged_bam": str(bam.resolve())}),
                encoding="utf-8",
            )

            sv = scratch_directory / "results" / "sv"
            sv.mkdir(parents=True)
            for suffix in (
                ".sawfish.vcf.gz",
                ".sawfish.vcf.gz.tbi",
                ".supporting_reads.json.gz",
            ):
                (sv / f"sample{suffix}").write_text("sv\n", encoding="utf-8")
            svtopo = sv / "sample.svtopo"
            svtopo.mkdir()
            (svtopo / "summary.tsv").write_text("svtopo\n", encoding="utf-8")

            hla = scratch_directory / "results" / "hla"
            hla.mkdir(parents=True)
            for suffix in (
                "_hifihla_report.json",
                "_hifihla_report.tsv",
                "_hifihla_summary.tsv",
            ):
                (hla / f"sample{suffix}").write_text("hla\n", encoding="utf-8")

            target = (
                run_directory
                / "run-metadata"
                / "sample-publication"
                / "sample.json"
            ).resolve()
            completed = subprocess.run(
                [
                    "snakemake",
                    str(target),
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
            published = run_directory / "samples" / "sample" / "pbWGS"
            manifest = json.loads((published / "outputs.json").read_text())
            self.assertEqual(
                Path(manifest["haplotagged_bam"]),
                (published / "sample.GRCh38.haplotagged.bam").resolve(),
            )
            self.assertEqual(json.loads(target.read_text())["status"], "complete")

    def test_fiber_branch_starts_with_input_detection(self):
        with tempfile.TemporaryDirectory() as temporary:
            completed = self._dry_run(
                Path(temporary),
                "auto-run",
                "yes",
                fiberseq_input="auto",
            )

        self.assertIn("checkpoint detect_fiberseq_input", completed.stdout)
        snakefile = (REPOSITORY_ROOT / "Snakefile").read_text(encoding="utf-8")
        predict_rule = snakefile.split("rule predict_m6a:", 1)[1].split(
            "\n\nrule ", 1
        )[0]
        self.assertIn("{params.executable:q} predict-m6a", predict_rule)
        self.assertIn('"envs/fibertools.yaml"', predict_rule)
        self.assertIn("{params.executable:q} --version", predict_rule)
        self.assertIn(
            'gpu=1 if FT_PREDICT_BACKEND == "external" else 0',
            predict_rule,
        )

    def test_legacy_fiber_runsheet_defaults_to_auto_detection(self):
        with tempfile.TemporaryDirectory() as temporary:
            completed = self._dry_run(
                Path(temporary),
                "legacy-fiber-run",
                "yes",
                fiberseq_input=None,
            )

        self.assertIn("checkpoint detect_fiberseq_input", completed.stdout)

    def test_detection_report_selects_only_the_required_preparation(self):
        cases = {
            "predict_m6a": "predict_m6a",
            "add_nucleosomes": "add_nucleosomes",
            "use_existing": None,
        }
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for action, expected_rule in cases.items():
                with self.subTest(action=action):
                    completed = self._dry_run(
                        root,
                        f"detected-{action}",
                        "yes",
                        detected_action=action,
                    )
                    for rule in ("predict_m6a", "add_nucleosomes"):
                        if rule == expected_rule:
                            self.assertIn(f"rule {rule}", completed.stdout)
                        else:
                            self.assertNotIn(f"rule {rule}", completed.stdout)

    def test_aligned_recovery_bam_runs_validation_and_current_wdl(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (
                run_directory,
                scratch_directory,
                _,
                config_path,
                environment,
            ) = self._fixture(
                root,
                "aligned-recovery",
                "yes",
                input_stage="aligned",
            )
            validation_run = self._run_dry_run(
                run_directory,
                config_path,
                environment,
                printshellcmds=True,
            )
            self._write_checkpoint_reports(scratch_directory)

            completed = self._run_dry_run(
                run_directory,
                config_path,
                environment,
            )

        self.assertIn("rule validate_reanalysis_input", validation_run.stdout)
        self.assertRegex(validation_run.stdout, r"(?m)^\s*outputs_json=$")
        self.assertIn('if [[ -n "$outputs_json" ]]', validation_run.stdout)
        self.assertIn("run_validator\n", validation_run.stdout)
        self.assertIn("rule hifi_wgs", completed.stdout)

    def test_reused_pbwgs_outputs_skip_current_wdl(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (
                run_directory,
                scratch_directory,
                _,
                config_path,
                environment,
            ) = self._fixture(
                root,
                "pbwgs-reuse",
                "yes",
                input_stage="pbwgs",
            )
            validation_run = self._run_dry_run(
                run_directory,
                config_path,
                environment,
                printshellcmds=True,
            )
            self._write_checkpoint_reports(scratch_directory)

            completed = self._run_dry_run(
                run_directory,
                config_path,
                environment,
            )

        self.assertIn("rule validate_reanalysis_input", validation_run.stdout)
        self.assertIn(
            "outputs_json="
            f"{(run_directory / 'inputs' / 'existing-outputs.json').resolve()}",
            validation_run.stdout,
        )
        self.assertIn('if [[ -n "$outputs_json" ]]', validation_run.stdout)
        self.assertIn(
            'run_validator --outputs-json "$outputs_json"',
            validation_run.stdout,
        )
        self.assertIn("rule copy_pbwdl", completed.stdout)
        self.assertNotIn("rule hifi_wgs", completed.stdout)

    def test_reanalysis_outputs_json_shell_handles_optional_path(self):
        snakefile = (REPOSITORY_ROOT / "Snakefile").read_text(encoding="utf-8")
        validation_rule = snakefile.split(
            "rule validate_reanalysis_input:", 1
        )[1].split("\n\nrule ", 1)[0]
        self.assertIn("outputs_json={params.outputs_json:q}", validation_rule)
        self.assertIn("run_validator()", validation_rule)
        self.assertIn('if [[ -n "$outputs_json" ]]', validation_rule)
        self.assertIn(
            'run_validator --outputs-json "$outputs_json"',
            validation_rule,
        )

        shell = r'''set -euo pipefail
outputs_json=$1
run_validator() {
    for argument in "$@"; do
        printf '%s\0' "$argument"
    done
}
if [[ -n "$outputs_json" ]]; then
    run_validator --outputs-json "$outputs_json"
else
    run_validator
fi
'''
        cases = {
            "aligned": ("", b""),
            "pbwgs": (
                "/tmp/path with spaces/outputs.json",
                b"--outputs-json\0/tmp/path with spaces/outputs.json\0",
            ),
        }
        for stage, (outputs_json, expected) in cases.items():
            with self.subTest(stage=stage):
                completed = subprocess.run(
                    ["/bin/bash", "-c", shell, "bash", outputs_json],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    check=False,
                )
                self.assertEqual(completed.returncode, 0, completed.stderr)
                self.assertEqual(completed.stdout, expected)

    def test_fire_only_reuses_pbwgs_bam_without_copying_all_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (
                run_directory,
                scratch_directory,
                _,
                config_path,
                environment,
            ) = self._fixture(
                root,
                "pbwgs-fire-reuse",
                "yes",
                input_stage="pbwgs",
            )
            self._write_checkpoint_reports(scratch_directory)
            self._write_eligibility(scratch_directory)

            completed = self._run_dry_run(
                run_directory,
                config_path,
                environment,
                target="fire_only",
            )

        self.assertIn("rule fire_input_bam", completed.stdout)
        self.assertNotIn("rule copy_pbwdl", completed.stdout)
        self.assertNotIn("rule hifi_wgs", completed.stdout)

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

            detection = (
                scratch_directory
                / "intermediate"
                / "fiberseq"
                / "input-detection"
                / "sample.json"
            )
            detection.parent.mkdir(parents=True)
            detection.write_text(
                json.dumps({"selected_action": "use_existing"}),
                encoding="utf-8",
            )

            self._write_eligibility(scratch_directory)

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
            for rule in (
                "fire_shuffled_pileup_chromosome",
                "fire_pileup_chromosome",
            ):
                rule_output = completed.stdout.split(f"rule {rule}:", 1)[1]
                rule_output = rule_output.split("\n[", 1)[0]
                with self.subTest(indexed_rule=rule):
                    self.assertIn("filtered.cram.crai", rule_output)
            self.assertNotIn("pixi", completed.stdout.lower())
            self.assertNotIn("rule scatac_map", completed.stdout)
            self.assertNotIn("rule sawfish", completed.stdout)
            self.assertNotIn("rule hifihla", completed.stdout)

    def _fixture(
        self,
        root: Path,
        run_name: str,
        hia5: str,
        fiberseq_input: str | None = "auto",
        input_stage: str = "raw",
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

        pbwgs_outputs = inputs / "existing-outputs.json"
        pbwgs_outputs.write_text(
            json.dumps(
                {
                    "humanwgs_singleton.workflow_version": "v2.0.7",
                    "humanwgs_singleton.stat_mean_depth": 20.0,
                    "humanwgs_singleton.merged_haplotagged_bam": str(
                        (inputs / "sample.GRCh38.haplotagged.bam").resolve()
                    ),
                }
            ),
            encoding="utf-8",
        )
        (inputs / "sample.GRCh38.haplotagged.bam").touch()

        runsheet = invocation / "runsheet.csv"
        if fiberseq_input is None:
            runsheet_text = (
                "core_name,lab_name,run,cell,barcode,hifi_bam,sex,hia5\n"
                f"sample,,,,,inputs/sample.bam,F,{hia5}\n"
            )
        else:
            hifi_bam = "" if input_stage == "pbwgs" else "inputs/sample.bam"
            outputs_json = (
                "inputs/existing-outputs.json"
                if input_stage == "pbwgs"
                else ""
            )
            runsheet_text = (
                "core_name,lab_name,run,cell,barcode,hifi_bam,input_stage,"
                "pbwgs_outputs_json,sex,hia5,fiberseq_input\n"
                f"sample,,,,,{hifi_bam},{input_stage},{outputs_json},F,"
                f"{hia5},{fiberseq_input}\n"
            )
        runsheet.write_text(runsheet_text, encoding="utf-8")
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
            "ft_predict_executable": "ft",
            "fiberseq_detection_reads": 5000,
            "fiberseq_detection_min_fraction": 0.5,
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

    def _write_checkpoint_reports(self, scratch_directory: Path) -> None:
        validation = (
            scratch_directory
            / "intermediate"
            / "reanalysis"
            / "input-validation"
            / "sample.json"
        )
        validation.parent.mkdir(parents=True)
        validation.write_text(
            json.dumps({"validation": "passed"}),
            encoding="utf-8",
        )
        detection = (
            scratch_directory
            / "intermediate"
            / "fiberseq"
            / "input-detection"
            / "sample.json"
        )
        detection.parent.mkdir(parents=True)
        detection.write_text(
            json.dumps({"selected_action": "use_existing"}),
            encoding="utf-8",
        )

    def _write_eligibility(self, scratch_directory: Path) -> None:
        manifest = (
            scratch_directory
            / "intermediate"
            / "fire"
            / "eligibility"
            / "sample.json"
        )
        manifest.parent.mkdir(parents=True)
        manifest.write_text(
            json.dumps(
                {
                    "sample": "sample",
                    "mean_depth": 20.0,
                    "minimum_mean_depth": 8.0,
                    "eligible": True,
                    "reason": "eligible",
                }
            ),
            encoding="utf-8",
        )
        report = manifest.with_suffix(".tsv")
        report.parent.mkdir(parents=True, exist_ok=True)
        report.write_text(
            "sample\tmean_depth\tminimum_mean_depth\teligible\treason\n"
            "sample\t20\t8\tyes\teligible\n",
            encoding="utf-8",
        )

    def _run_dry_run(
        self,
        run_directory: Path,
        config_path: Path,
        environment: dict[str, str],
        target: str | None = None,
        printshellcmds: bool = False,
    ) -> subprocess.CompletedProcess[str]:
        command = ["snakemake"]
        if target is not None:
            command.append(target)
        command.extend(
            (
                "--snakefile",
                str(REPOSITORY_ROOT / "Snakefile"),
                "--configfile",
                str(config_path),
                "--dry-run",
                "--cores",
                "1",
            )
        )
        if printshellcmds:
            command.append("--printshellcmds")
        completed = subprocess.run(
            command,
            cwd=run_directory,
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        self.assertEqual(completed.returncode, 0, completed.stdout)
        return completed

    def _dry_run(
        self,
        root: Path,
        run_name: str,
        hia5: str,
        fiberseq_input: str | None = "auto",
        detected_action: str | None = None,
    ) -> subprocess.CompletedProcess[str]:
        (
            run_directory,
            scratch_directory,
            _,
            config_path,
            environment,
        ) = self._fixture(root, run_name, hia5, fiberseq_input)
        if detected_action is not None:
            detection = (
                scratch_directory
                / "intermediate"
                / "fiberseq"
                / "input-detection"
                / "sample.json"
            )
            detection.parent.mkdir(parents=True)
            detection.write_text(
                json.dumps({"selected_action": detected_action}),
                encoding="utf-8",
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
            if detected_action is None:
                self.assertIn(
                    "checkpoint detect_fiberseq_input",
                    completed.stdout,
                )
            else:
                self.assertNotIn("checkpoint detect_fiberseq_input", completed.stdout)
        else:
            self.assertNotIn("checkpoint fire_eligibility", completed.stdout)
            self.assertNotIn("checkpoint detect_fiberseq_input", completed.stdout)
            self.assertNotIn("rule add_nucleosomes", completed.stdout)
            self.assertNotIn("rule predict_m6a", completed.stdout)
        return completed
