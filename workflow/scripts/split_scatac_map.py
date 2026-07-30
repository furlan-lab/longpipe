from contextlib import ExitStack
from pathlib import Path


output_paths = [
    Path(str(snakemake.output.maximum)),
    Path(str(snakemake.output.median)),
    Path(str(snakemake.output.total)),
]
for output_path in output_paths:
    output_path.parent.mkdir(parents=True, exist_ok=True)

record_count = 0
with ExitStack() as stack:
    input_handle = stack.enter_context(
        open(str(snakemake.input.mapped), encoding="utf-8")
    )
    output_handles = [
        stack.enter_context(open(output_path, "w", encoding="utf-8"))
        for output_path in output_paths
    ]

    for line_number, line in enumerate(input_handle, start=1):
        stripped_line = line.rstrip("\n")
        if not stripped_line:
            continue
        if stripped_line.startswith("#"):
            for output_handle in output_handles:
                output_handle.write(line)
            continue

        fields = stripped_line.split("\t")
        if len(fields) < 6:
            raise ValueError(
                f"{snakemake.input.mapped}:{line_number}: "
                "expected BED columns followed by max, median, and sum"
            )

        interval_fields = fields[:-3]
        metrics = fields[-3:]
        for output_handle, metric in zip(output_handles, metrics):
            output_handle.write("\t".join([*interval_fields, metric]) + "\n")
        record_count += 1

if record_count == 0:
    raise ValueError(f"{snakemake.input.mapped}: no interval records found")
