import argparse
import gzip
import hashlib
import itertools
from contextlib import contextmanager
from pathlib import Path


@contextmanager
def open_decompressed(path):
    with path.open("rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    with opener(path, "rb") as handle:
        yield handle


def decompressed_digest(path):
    digest = hashlib.sha256()
    with open_decompressed(path) as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def compare_pileups(first_path, second_path):
    with open_decompressed(first_path) as first_handle:
        with open_decompressed(second_path) as second_handle:
            lines = itertools.zip_longest(first_handle, second_handle)
            for line_number, (first_line, second_line) in enumerate(
                lines, start=1
            ):
                if first_line != second_line:
                    return line_number
    return None


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("legacy_pileup", type=Path)
    parser.add_argument("native_pileup", type=Path)
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    differing_line = compare_pileups(
        arguments.legacy_pileup,
        arguments.native_pileup,
    )
    if differing_line is not None:
        print(f"pileups differ at decompressed line {differing_line}")
        raise SystemExit(1)
    digest = decompressed_digest(arguments.legacy_pileup)
    print(f"pileups are identical after decompression: sha256={digest}")


if __name__ == "__main__":
    main()
