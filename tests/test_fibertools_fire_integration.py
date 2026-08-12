from __future__ import annotations

import os
import subprocess
import tempfile
import unittest
from pathlib import Path


FT = os.environ.get("LONGPIPE_TEST_FT")
LEGACY_FIXTURE = os.environ.get("LONGPIPE_TEST_FIBERTOOLS_BAM")
SAMTOOLS = os.environ.get("LONGPIPE_TEST_SAMTOOLS", "samtools")


@unittest.skipUnless(
    FT and LEGACY_FIXTURE,
    "set LONGPIPE_TEST_FT and LONGPIPE_TEST_FIBERTOOLS_BAM for integration",
)
class FibertoolsFireIntegrationTest(unittest.TestCase):
    def test_fire_reads_ma_nucleosome_and_msp_annotations(self):
        assert FT is not None
        assert LEGACY_FIXTURE is not None
        fixture = Path(LEGACY_FIXTURE)
        self.assertTrue(fixture.is_file(), fixture)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ma_bam = root / "ma.bam"
            annotations_path = root / "annotations.tsv"
            features = root / "fire-features.tsv"

            converted = subprocess.run(
                [FT, "convert-tags", str(fixture), str(ma_bam)],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False,
            )
            self.assertEqual(converted.returncode, 0, converted.stdout)

            converted_sam = subprocess.run(
                [SAMTOOLS, "view", str(ma_bam)],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False,
            )
            self.assertEqual(converted_sam.returncode, 0, converted_sam.stdout)
            records = converted_sam.stdout.splitlines()
            self.assertTrue(records, converted_sam.stdout)
            self.assertTrue(
                any("\tMA:Z:" in f"\t{record}" for record in records),
                "converted BAM has no MA annotations",
            )
            for legacy_tag in ("ns", "nl", "as", "al"):
                self.assertFalse(
                    any(f"\t{legacy_tag}:" in f"\t{record}" for record in records),
                    f"legacy {legacy_tag} annotation survived conversion",
                )

            annotations = subprocess.run(
                [
                    FT,
                    "extract",
                    "--all",
                    str(annotations_path),
                    str(ma_bam),
                ],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False,
            )
            self.assertEqual(annotations.returncode, 0, annotations.stdout)
            header, *rows = annotations_path.read_text(
                encoding="utf-8"
            ).splitlines()
            self.assertIn("nuc_starts", header)
            self.assertIn("msp_starts", header)
            self.assertTrue(rows, annotations.stdout)

            fire = subprocess.run(
                [
                    FT,
                    "fire",
                    "--feats-to-text",
                    str(ma_bam),
                    str(features),
                ],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False,
            )
            self.assertEqual(fire.returncode, 0, fire.stdout)
            feature_lines = features.read_text(encoding="utf-8").splitlines()
            self.assertGreater(len(feature_lines), 1, fire.stdout)
            self.assertIn("msp_len", feature_lines[0].split("\t"))


if __name__ == "__main__":
    unittest.main()
