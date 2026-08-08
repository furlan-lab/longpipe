#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate an aligned BAM used for Longpipe reanalysis."
    )
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--reference-index", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--stage", required=True, choices=("aligned", "pbwgs"))
    parser.add_argument("--outputs-json", type=Path)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--records", type=int, default=5000)
    return parser.parse_args()


def fail(message: str) -> None:
    raise ValueError(message)


def run_checked(command: list[str]) -> str:
    completed = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if completed.returncode != 0:
        detail = completed.stderr.strip() or completed.stdout.strip()
        raise ValueError(f"command failed ({' '.join(command)}): {detail}")
    return completed.stdout


def parse_bam_header(header: str) -> tuple[str | None, dict[str, int]]:
    sort_order = None
    contigs: dict[str, int] = {}
    for line in header.splitlines():
        fields = line.split("\t")
        if fields[0] == "@HD":
            for field in fields[1:]:
                if field.startswith("SO:"):
                    sort_order = field[3:]
        elif fields[0] == "@SQ":
            values = {
                field[:2]: field[3:]
                for field in fields[1:]
                if len(field) > 3 and field[2] == ":"
            }
            try:
                name = values["SN"]
                length = int(values["LN"])
            except (KeyError, ValueError) as error:
                raise ValueError(f"malformed @SQ header line: {line}") from error
            if name in contigs:
                raise ValueError(f"duplicate BAM reference contig: {name}")
            contigs[name] = length
    return sort_order, contigs


def read_reference_index(path: Path) -> dict[str, int]:
    contigs: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                fail(f"{path}:{line_number}: malformed FASTA index line")
            try:
                length = int(fields[1])
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: invalid contig length"
                ) from error
            if fields[0] in contigs:
                fail(f"{path}:{line_number}: duplicate contig {fields[0]}")
            contigs[fields[0]] = length
    if not contigs:
        fail(f"reference index contains no contigs: {path}")
    return contigs


def sample_records(
    samtools: str, bam: Path, maximum: int
) -> dict[str, int]:
    if maximum < 1:
        fail("--records must be at least one")
    process = subprocess.Popen(
        [samtools, "view", str(bam)],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert process.stdout is not None
    counts = {
        "records": 0,
        "primary": 0,
        "mapped_primary": 0,
        "secondary": 0,
        "supplementary": 0,
    }
    try:
        for line in process.stdout:
            fields = line.split("\t", 3)
            if len(fields) < 2:
                fail("samtools returned a malformed alignment record")
            try:
                flag = int(fields[1])
            except ValueError as error:
                raise ValueError("samtools returned an invalid BAM flag") from error
            counts["records"] += 1
            if flag & 0x100:
                counts["secondary"] += 1
            elif flag & 0x800:
                counts["supplementary"] += 1
            else:
                counts["primary"] += 1
                if not flag & 0x4:
                    counts["mapped_primary"] += 1
            if counts["records"] >= maximum:
                break
    finally:
        if process.poll() is None:
            process.terminate()
        process.communicate()
    if counts["records"] == 0:
        fail(f"BAM contains no alignment records: {bam}")
    if counts["mapped_primary"] == 0:
        fail(
            f"sampled {counts['records']} BAM records but found no mapped "
            "primary reads"
        )
    return counts


def workflow_version(path: Path | None) -> object:
    if path is None:
        return None
    with path.open(encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        fail(f"PacBio WDL outputs are not a JSON object: {path}")
    outputs = payload.get("outputs", payload)
    if not isinstance(outputs, dict):
        fail(f"PacBio WDL outputs are not a JSON object: {path}")
    return outputs.get("humanwgs_singleton.workflow_version")


def atomic_write(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)
            handle.write("\n")
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def main() -> int:
    arguments = parse_arguments()
    if not arguments.input.is_file():
        fail(f"input BAM does not exist: {arguments.input}")
    if not arguments.reference_index.is_file():
        fail(f"reference index does not exist: {arguments.reference_index}")
    if arguments.stage == "pbwgs" and arguments.outputs_json is None:
        fail("input_stage=pbwgs requires --outputs-json")
    if arguments.stage == "aligned" and arguments.outputs_json is not None:
        fail("input_stage=aligned cannot use --outputs-json")

    quickcheck = subprocess.run(
        [arguments.samtools, "quickcheck", "-v", str(arguments.input)],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if quickcheck.returncode != 0:
        detail = quickcheck.stdout.strip() or quickcheck.stderr.strip()
        fail(f"samtools quickcheck failed for {arguments.input}: {detail}")

    header = run_checked(
        [arguments.samtools, "view", "-H", str(arguments.input)]
    )
    sort_order, bam_contigs = parse_bam_header(header)
    reference_contigs = read_reference_index(arguments.reference_index)
    if sort_order != "coordinate":
        fail(
            f"{arguments.sample}: reanalysis BAM must be coordinate sorted; "
            f"@HD SO is {sort_order!r}"
        )
    bam_dictionary = list(bam_contigs.items())
    reference_dictionary = list(reference_contigs.items())
    if bam_dictionary != reference_dictionary:
        missing = sorted(set(reference_contigs) - set(bam_contigs))
        unexpected = sorted(set(bam_contigs) - set(reference_contigs))
        length_mismatches = sorted(
            name
            for name in set(bam_contigs) & set(reference_contigs)
            if bam_contigs[name] != reference_contigs[name]
        )
        order_mismatch = (
            not missing
            and not unexpected
            and not length_mismatches
            and list(bam_contigs) != list(reference_contigs)
        )
        fail(
            f"{arguments.sample}: BAM reference dictionary does not match "
            f"{arguments.reference_index}; missing={missing[:5]}, "
            f"unexpected={unexpected[:5]}, "
            f"length_mismatches={length_mismatches[:5]}, "
            f"order_mismatch={order_mismatch}"
        )

    sampled = sample_records(
        arguments.samtools,
        arguments.input,
        arguments.records,
    )
    report = {
        "sample": arguments.sample,
        "input_stage": arguments.stage,
        "input_bam": str(arguments.input.resolve()),
        "outputs_json": (
            str(arguments.outputs_json.resolve())
            if arguments.outputs_json is not None
            else None
        ),
        "source_workflow_version": workflow_version(arguments.outputs_json),
        "reference_index": str(arguments.reference_index.resolve()),
        "reference_contigs": len(reference_contigs),
        "sort_order": sort_order,
        "sampled_records": sampled,
        "reference_dictionary_matches": True,
        "validation": "passed",
        "action": (
            "reuse_pbwgs_outputs"
            if arguments.stage == "pbwgs"
            else "realign_with_current_pbwgs"
        ),
    }
    atomic_write(arguments.output, report)
    print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
