from __future__ import annotations

import csv
import gzip
import importlib.util
import tempfile
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def load_script(name: str):
    path = REPOSITORY_ROOT / "workflow" / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


haplotype_tracks = load_script("fire_haplotype_tracks")
haplotype_differences = load_script("fire_haplotype_differences")


class FireHaplotypeTracksTest(unittest.TestCase):
    def test_tracks_report_percent_accessible_and_haplotype_difference(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            pileup = root / "pileup.bed.gz"
            with gzip.open(pileup, "wt", encoding="utf-8", newline="") as handle:
                handle.write(
                    "#chrom\tstart\tend\tcoverage_H1\tfire_coverage_H1\t"
                    "coverage_H2\tfire_coverage_H2\n"
                    "chr1\t0\t10\t10\t4\t10\t6\n"
                    "chr1\t10\t20\t0\t0\t5\t1\n"
                    "chr1\t20\t30\t.\t.\t4\t2\n"
                )
            h1 = root / "H1.bedgraph"
            h2 = root / "H2.bedgraph"
            difference = root / "difference.bedgraph"

            haplotype_tracks.build_tracks(pileup, h1, h2, difference)

            self.assertEqual(h1.read_text(), "chr1\t0\t10\t40\n")
            self.assertEqual(
                h2.read_text(),
                "chr1\t0\t10\t60\n"
                "chr1\t10\t20\t20\n"
                "chr1\t20\t30\t50\n",
            )
            self.assertEqual(difference.read_text(), "chr1\t0\t10\t-20\n")


class FireHaplotypeDifferencesTest(unittest.TestCase):
    def test_fisher_exact_matches_reference_value(self):
        pvalue = haplotype_differences.fisher_exact_two_sided(1, 9, 11, 3)

        self.assertAlmostEqual(pvalue, 0.0027594561852200836, places=14)

    def test_benjamini_hochberg_preserves_input_order(self):
        adjusted = haplotype_differences.benjamini_hochberg(
            [0.01, 0.04, 0.03]
        )

        self.assertEqual(adjusted, [0.03, 0.04, 0.04])

    def test_analysis_filters_and_summarizes_confident_local_maxima(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            pileup = root / "pileup.bed"
            pileup.write_text(
                "#chrom\tstart\tend\tFDR\tis_local_max\tcoverage_H1\t"
                "fire_coverage_H1\tcoverage_H2\tfire_coverage_H2\n"
                "chr1\t0\t1\t0.01\ttrue\t10\t1\t14\t11\n"
                "chr1\t10\t11\t0.02\ttrue\t10\t5\t10\t5\n"
                "chr1\t20\t21\t0.01\tfalse\t10\t1\t14\t11\n"
                "chr1\t30\t31\t0.10\ttrue\t10\t1\t14\t11\n"
                "chr1\t40\t41\t0.01\ttrue\t3\t1\t14\t11\n",
                encoding="utf-8",
            )
            output = root / "differences.bed"
            summary = root / "summary.tsv"

            haplotype_differences.analyze_pileup(
                pileup,
                output,
                summary,
                minimum_coverage=5,
                maximum_fire_fdr=0.05,
                qvalue_threshold=0.05,
            )

            with output.open(newline="", encoding="utf-8") as handle:
                records = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(records), 2)
            self.assertEqual(records[0]["#chrom"], "chr1")
            self.assertEqual(records[0]["significant"], "yes")
            self.assertAlmostEqual(float(records[0]["p_value"]), 0.002759456185)
            self.assertAlmostEqual(float(records[0]["q_value"]), 0.00551891237)
            self.assertEqual(records[1]["significant"], "no")

            with summary.open(newline="", encoding="utf-8") as handle:
                metrics = dict(csv.reader(handle, delimiter="\t"))
            self.assertEqual(metrics["candidate_local_maxima"], "3")
            self.assertEqual(metrics["tested_loci"], "2")
            self.assertEqual(metrics["significant_loci"], "1")


if __name__ == "__main__":
    unittest.main()
