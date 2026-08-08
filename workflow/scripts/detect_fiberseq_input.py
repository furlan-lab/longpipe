#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path


KINETICS_TAGS = ("fi", "ri", "fp", "rp")
LEGACY_ANNOTATION_TAGS = ("ns", "nl", "as", "al")
MA_ANNOTATION_TYPES = ("nuc", "msp")
M6A_PREFIX = re.compile(r"A[+-]a[.?]?")


class DetectionError(RuntimeError):
    pass


def parse_tags(fields: list[str]) -> dict[str, tuple[str, str]]:
    tags: dict[str, tuple[str, str]] = {}
    for field in fields:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = (parts[1], parts[2])
    return tags


def array_tag_present(tags: dict[str, tuple[str, str]], name: str) -> bool:
    value = tags.get(name)
    return value is not None and value[0] == "B" and bool(value[1])


def has_m6a_calls(mm_value: str) -> bool:
    for group in mm_value.rstrip(";").split(";"):
        prefix, separator, calls = group.partition(",")
        if separator and M6A_PREFIX.fullmatch(prefix) and re.search(r"\d", calls):
            return True
    return False


def has_ma_footprints(ma_value: str) -> bool:
    annotation_types: set[str] = set()
    for section in ma_value.split(";")[1:]:
        type_info, separator, _ = section.partition(":")
        if not separator:
            continue
        for annotation_type in MA_ANNOTATION_TYPES:
            if type_info.startswith(
                (
                    f"{annotation_type}+",
                    f"{annotation_type}-",
                    f"{annotation_type}.",
                )
            ):
                annotation_types.add(annotation_type)
    return all(name in annotation_types for name in MA_ANNOTATION_TYPES)


def inspect_records(
    bam: Path,
    max_reads: int,
    samtools: str,
) -> dict[str, int]:
    command = [samtools, "view", "-F", "2304", str(bam)]
    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except OSError as error:
        raise DetectionError(f"could not run {samtools!r}: {error}") from error

    counts = {
        "reads_tested": 0,
        "m6a_mm_reads": 0,
        "m6a_reads": 0,
        "m6a_without_ml_reads": 0,
        "annotation_reads": 0,
        "ma_annotation_reads": 0,
        "legacy_annotation_reads": 0,
        "kinetics_reads": 0,
    }
    stopped_early = False
    try:
        assert process.stdout is not None
        for line in process.stdout:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                raise DetectionError(
                    f"samtools returned a malformed alignment for {bam}"
                )
            tags = parse_tags(fields[11:])
            mm = tags.get("MM") or tags.get("Mm")
            ml = tags.get("ML") or tags.get("Ml")
            has_m6a = mm is not None and has_m6a_calls(mm[1])

            counts["reads_tested"] += 1
            if has_m6a:
                counts["m6a_mm_reads"] += 1
                if ml is not None and array_tag_present(
                    {"ML": ml}, "ML"
                ):
                    counts["m6a_reads"] += 1
                else:
                    counts["m6a_without_ml_reads"] += 1
            ma = tags.get("MA")
            has_ma_annotations = (
                ma is not None
                and ma[0] == "Z"
                and has_ma_footprints(ma[1])
            )
            has_legacy_annotations = all(
                array_tag_present(tags, tag)
                for tag in LEGACY_ANNOTATION_TAGS
            )
            if has_ma_annotations:
                counts["ma_annotation_reads"] += 1
            if has_legacy_annotations:
                counts["legacy_annotation_reads"] += 1
            if has_ma_annotations or has_legacy_annotations:
                counts["annotation_reads"] += 1
            if all(array_tag_present(tags, tag) for tag in KINETICS_TAGS):
                counts["kinetics_reads"] += 1

            if counts["reads_tested"] >= max_reads:
                stopped_early = True
                process.terminate()
                break

        _, stderr = process.communicate()
        if not stopped_early and process.returncode:
            detail = stderr.strip() or f"exit status {process.returncode}"
            raise DetectionError(f"samtools could not inspect {bam}: {detail}")
    finally:
        if process.poll() is None:
            process.kill()
            process.wait()

    if counts["reads_tested"] == 0:
        raise DetectionError(f"no primary alignments found in {bam}")
    return counts


def select_action(
    counts: dict[str, int],
    mode: str,
    min_fraction: float,
) -> tuple[str, str, dict[str, float]]:
    total = counts["reads_tested"]
    fractions = {
        "m6a_fraction": counts["m6a_reads"] / total,
        "annotation_fraction": counts["annotation_reads"] / total,
        "kinetics_fraction": counts["kinetics_reads"] / total,
    }

    has_m6a = fractions["m6a_fraction"] >= min_fraction
    has_annotations = fractions["annotation_fraction"] >= min_fraction
    has_kinetics = fractions["kinetics_fraction"] >= min_fraction

    if mode == "kinetics":
        if not has_kinetics:
            raise DetectionError(
                "fiberseq_input=kinetics requires fi, ri, fp, and rp on at "
                f"least {min_fraction:.0%} of reads; observed "
                f"{fractions['kinetics_fraction']:.1%}"
            )
        return "predict_m6a", "explicit_kinetics_override", fractions

    if mode == "m6a":
        if not has_m6a:
            raise DetectionError(
                "fiberseq_input=m6a requires A+a MM/ML calls on at least "
                f"{min_fraction:.0%} of reads; observed "
                f"{fractions['m6a_fraction']:.1%}"
            )
        if has_annotations:
            return "use_existing", "explicit_m6a_already_annotated", fractions
        return "add_nucleosomes", "explicit_m6a_override", fractions

    if has_m6a:
        if has_annotations:
            return "use_existing", "m6a_and_footprints_detected", fractions
        return "add_nucleosomes", "m6a_detected_without_footprints", fractions

    if counts["m6a_mm_reads"]:
        raise DetectionError(
            "m6A calls were present on too few reads for automatic selection: "
            f"{fractions['m6a_fraction']:.1%}; set an explicit override only "
            "after checking the BAM"
        )
    if has_kinetics:
        return "predict_m6a", "kinetics_detected_without_m6a", fractions

    raise DetectionError(
        "BAM has neither sufficient A+a MM/ML calls nor sufficient PacBio "
        "fi/ri/fp/rp kinetics for Fiber-seq processing "
        f"(m6A {fractions['m6a_fraction']:.1%}, kinetics "
        f"{fractions['kinetics_fraction']:.1%})"
    )


def write_report(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
        json.dump(payload, handle, indent=2)
        handle.write("\n")
    try:
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Detect the required Fiber-seq BAM preparation step"
    )
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument(
        "--mode",
        choices=("auto", "m6a", "kinetics"),
        default="auto",
    )
    parser.add_argument("--reads", type=int, default=5000)
    parser.add_argument("--min-fraction", type=float, default=0.5)
    parser.add_argument("--samtools", default="samtools")
    args = parser.parse_args()
    if args.reads < 1:
        parser.error("--reads must be at least 1")
    if not 0 < args.min_fraction <= 1:
        parser.error("--min-fraction must be greater than 0 and at most 1")
    return args


def main() -> int:
    args = parse_args()
    try:
        counts = inspect_records(args.input, args.reads, args.samtools)
        action, reason, fractions = select_action(
            counts,
            args.mode,
            args.min_fraction,
        )
        payload: dict[str, object] = {
            "sample": args.sample,
            "input_bam": str(args.input.resolve()),
            "requested_mode": args.mode,
            "selected_action": action,
            "decision_reason": reason,
            "maximum_reads": args.reads,
            "minimum_fraction": args.min_fraction,
            **counts,
            **fractions,
            "kinetics_tags": list(KINETICS_TAGS),
            "ma_annotation_tag": "MA",
            "ma_annotation_types": list(MA_ANNOTATION_TYPES),
            "legacy_annotation_tags": list(LEGACY_ANNOTATION_TAGS),
        }
        write_report(args.output, payload)
    except DetectionError as error:
        print(f"{args.sample}: {error}", file=sys.stderr)
        return 2

    print(
        f"{args.sample}: selected {action} from {counts['reads_tested']} reads "
        f"(m6A {fractions['m6a_fraction']:.1%}, annotations "
        f"{fractions['annotation_fraction']:.1%}, kinetics "
        f"{fractions['kinetics_fraction']:.1%})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
