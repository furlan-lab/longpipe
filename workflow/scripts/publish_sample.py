"""Snakemake entry point for per-sample publication."""

import traceback
from pathlib import Path

from longpipe.publication import publish_fire_only, publish_sample


common = {
    "sample": str(snakemake.wildcards.sample),
    "scratch_directory": Path(str(snakemake.params.scratch_directory)),
    "run_directory": Path(str(snakemake.params.run_directory)),
    "status_path": Path(str(snakemake.output.status)),
    "fire_version": str(snakemake.params.fire_version),
    "invocation_id": str(snakemake.params.invocation_id),
}

log_path = Path(str(snakemake.log[0]))
log_path.parent.mkdir(parents=True, exist_ok=True)
try:
    decision_value = str(snakemake.params.fire_decision)
    decision = Path(decision_value) if decision_value else None
    if snakemake.params.mode == "fire_only":
        assert decision is not None
        publish_fire_only(
            **common,
            fire_decision_path=decision,
        )
    else:
        publish_sample(
            **common,
            peak_sets=list(snakemake.params.peak_sets),
            fire_decision_path=decision,
        )
except BaseException:
    with log_path.open("w", encoding="utf-8") as handle:
        traceback.print_exc(file=handle)
    raise
else:
    log_path.write_text(
        f"Published {common['sample']} ({snakemake.params.mode})\n",
        encoding="utf-8",
    )
