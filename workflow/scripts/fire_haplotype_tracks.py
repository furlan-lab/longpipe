import argparse
import csv
import gzip
import math
from contextlib import contextmanager
from pathlib import Path


@contextmanager
def open_text(path):
    with path.open("rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        yield handle


def parse_number(value):
    if value in {None, "", ".", "NA", "NaN", "nan"}:
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def percent_accessible(row, haplotype):
    coverage = parse_number(row[f"coverage_{haplotype}"])
    fire_coverage = parse_number(row[f"fire_coverage_{haplotype}"])
    if coverage is None or fire_coverage is None or coverage <= 0:
        return None
    return 100.0 * fire_coverage / coverage


def write_track(handle, chromosome, start, end, value):
    handle.write(f"{chromosome}\t{start}\t{end}\t{value:.10g}\n")


def build_tracks(input_path, h1_path, h2_path, difference_path):
    required = {
        "#chrom",
        "start",
        "end",
        "coverage_H1",
        "fire_coverage_H1",
        "coverage_H2",
        "fire_coverage_H2",
    }
    with open_text(input_path) as input_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError(
                "FIRE pileup is missing haplotype columns: "
                + ", ".join(missing)
            )
        with (
            h1_path.open("w", encoding="utf-8") as h1_handle,
            h2_path.open("w", encoding="utf-8") as h2_handle,
            difference_path.open("w", encoding="utf-8") as difference_handle,
        ):
            for row in reader:
                chromosome = row["#chrom"]
                start = row["start"]
                end = row["end"]
                h1_value = percent_accessible(row, "H1")
                h2_value = percent_accessible(row, "H2")
                if h1_value is not None:
                    write_track(h1_handle, chromosome, start, end, h1_value)
                if h2_value is not None:
                    write_track(h2_handle, chromosome, start, end, h2_value)
                if h1_value is not None and h2_value is not None:
                    write_track(
                        difference_handle,
                        chromosome,
                        start,
                        end,
                        h1_value - h2_value,
                    )


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("pileup", type=Path)
    parser.add_argument("h1_bedgraph", type=Path)
    parser.add_argument("h2_bedgraph", type=Path)
    parser.add_argument("difference_bedgraph", type=Path)
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    build_tracks(
        arguments.pileup,
        arguments.h1_bedgraph,
        arguments.h2_bedgraph,
        arguments.difference_bedgraph,
    )


if __name__ == "__main__":
    main()
