import argparse
import gzip
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
import polars.selectors as selectors


ROLLING_FIRE_SCORE_WINDOW_SIZE = 200


def is_gzipped(path):
    with path.open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def find_nearest(array, value):
    indices = np.searchsorted(array, value, side="left")
    indices[indices < 0] = 0
    indices[indices >= len(array)] = len(array) - 1
    return indices


def read_csv(*args, **kwargs):
    try:
        return pl.read_csv(*args, **kwargs)
    except pl.exceptions.NoDataError as error:
        print(
            "No data was found in the FIRE input file.",
            file=sys.stderr,
        )
        print(error, file=sys.stderr)
        raise


def read_pileup_file(infile, nrows):
    header = read_csv(infile, separator="\t", n_rows=1).columns
    opener = gzip.open if is_gzipped(infile) else open
    line_index = -1
    with opener(infile, "rt", encoding="utf-8") as handle:
        for line_index, _line in enumerate(handle):
            if line_index > 1:
                break
    if line_index < 2:
        return None

    schema_overrides = {}
    for column_index, column_name in enumerate(header, start=1):
        positional_name = f"column_{column_index}"
        if column_name in (
            "score",
            "score_H1",
            "score_H2",
            "score_shuffled",
        ):
            schema_overrides[positional_name] = pl.Float64
        elif column_name == "#chrom":
            schema_overrides[positional_name] = pl.Utf8
    logging.info("Header of the pileup file:\n%s", header)
    logging.info("Schema overrides:\n%s", schema_overrides)
    pileup = read_csv(
        infile,
        separator="\t",
        has_header=False,
        new_columns=header,
        comment_prefix="#",
        n_rows=nrows,
        infer_schema_length=100000,
        schema_overrides=schema_overrides,
    )
    logging.info("Done reading pileup file:\n%s", pileup)
    return pileup


def fdr_from_fire_scores(fire_scores):
    shuffled_counts = []
    real_counts = []
    thresholds = []
    current_real_count = 0.0
    current_shuffled_count = 0.0
    previous_score = -1.0
    first_score = True
    for score, is_real, base_pairs in fire_scores.iter_rows():
        if score != previous_score and current_real_count > 0 and not first_score:
            real_counts.append(current_real_count)
            shuffled_counts.append(current_shuffled_count)
            thresholds.append(previous_score)
        if score < 0.0:
            break
        if is_real:
            current_real_count += base_pairs
        else:
            current_shuffled_count += base_pairs
        previous_score = score
        first_score = False

    real_counts.append(1)
    shuffled_counts.append(1)
    thresholds.append(-1.0)
    shuffled_counts = np.array(shuffled_counts)
    real_counts = np.array(real_counts)
    thresholds = np.array(thresholds)
    false_discovery_rates = shuffled_counts / real_counts
    false_discovery_rates[false_discovery_rates > 1] = 1.0
    return thresholds, false_discovery_rates, shuffled_counts, real_counts


def fdr_table_from_scores(fire_scores):
    (
        thresholds,
        false_discovery_rates,
        shuffled_counts,
        real_counts,
    ) = fdr_from_fire_scores(fire_scores)
    results = pd.DataFrame(
        {
            "threshold": thresholds,
            "FDR": false_discovery_rates,
            "shuffled_bp": shuffled_counts,
            "real_bp": real_counts,
        }
    )
    results = results.groupby("FDR", sort=False).tail(1).reset_index(drop=True)
    results = (
        results.groupby("shuffled_bp", sort=False)
        .tail(1)
        .reset_index(drop=True)
    )
    results = (
        results.groupby("real_bp", sort=False)
        .tail(1)
        .reset_index(drop=True)
    )
    results["threshold"] = results["threshold"].round(2)
    results = (
        results.groupby("threshold", sort=False)
        .tail(1)
        .reset_index(drop=True)
    )
    results = results.sort_values("threshold")
    logging.info("FDR results\n%s", results)
    return results


def make_fdr_table(
    infile,
    outfile,
    nrows,
    max_coverage=None,
    min_coverage=None,
    max_fdr=0.05,
):
    pileup = read_pileup_file(infile, nrows)
    if pileup is None:
        raise ValueError("FIRE shuffled pileup contains too few records")
    if max_coverage is not None:
        pileup = pileup.filter(
            pl.col("coverage") <= max_coverage,
            pl.col("coverage_shuffled") <= max_coverage,
        )
    if min_coverage is not None:
        pileup = pileup.filter(
            pl.col("coverage") >= min_coverage,
            pl.col("coverage_shuffled") >= min_coverage,
        )

    fire_scores = (
        pileup.melt(
            value_vars=["score", "score_shuffled"],
            id_vars=["#chrom", "start", "end"],
            variable_name="type",
            value_name="score",
        )
        .with_columns(
            pl.when(pl.col("type") == "score")
            .then(True)
            .otherwise(False)
            .alias("is_real"),
            pl.col("end").sub(pl.col("start")).alias("bp"),
        )
        .group_by(["score", "is_real"])
        .agg(pl.sum("bp").alias("bp"))
        .sort("score", descending=True)
    )
    sums = fire_scores.group_by("is_real").agg(
        pl.sum("bp").alias("Mbp") / 1_000_000
    )
    logging.info("Number of Mbp in each category:\n%s", sums)
    logging.info("Done aggregating pileup file:\n%s", fire_scores)
    fdr_table = fdr_table_from_scores(fire_scores)
    fdr_table.to_csv(outfile, sep="\t", index=False)
    if fdr_table["FDR"].min() > max_fdr:
        raise ValueError(
            f"No FIRE score threshold has an FDR below {max_fdr}."
        )
    return fdr_table


def read_fdr_table(infile):
    fdr_table = read_csv(infile, separator="\t").to_pandas()
    logging.info("Read FDR table:\n%s", fdr_table)
    return fdr_table


def apply_fdr_table(infile, outfile, fdr_table, nrows):
    pileup = read_pileup_file(infile, nrows)
    if pileup is None:
        outfile.touch()
        return

    logging.info("Applying FDR table to pileup file:\n%s", pileup)
    rolling_max_score = (
        pileup.rolling(
            index_column="start",
            period=f"{ROLLING_FIRE_SCORE_WINDOW_SIZE}i",
            offset=f"-{ROLLING_FIRE_SCORE_WINDOW_SIZE // 2}i",
            closed="both",
            group_by="#chrom",
        )
        .agg(
            pl.max("score").alias("max_window_score").fill_null(-1.0),
        )
        .drop("#chrom")
    )
    pileup = (
        pileup.with_columns(rolling_max_score)
        .with_columns(
            pl.when(
                (pl.col("score") == pl.col("max_window_score"))
                & (pl.col("score") > 0.0)
            )
            .then(True)
            .otherwise(False)
            .alias("is_local_max"),
        )
        .with_columns(
            pl.col("is_local_max").rle_id().alias("local_max_group"),
        )
    )
    middle = (
        pileup.group_by(["#chrom", "local_max_group"])
        .agg(
            pl.col("end").sub(pl.col("start")).sum().alias("width"),
            pl.min("start").alias("first_start"),
        )
        .with_columns(
            pl.col("first_start").add(pl.col("width") // 2).alias("middle"),
        )
        .drop("width", "first_start")
    )
    pileup = (
        pileup.join(middle, on=["#chrom", "local_max_group"], how="left")
        .with_columns(
            (
                (pl.col("start") <= pl.col("middle"))
                & (pl.col("end") > pl.col("middle"))
            ).alias("is_middle"),
        )
        .with_columns(
            (pl.col("is_local_max") & pl.col("is_middle")).alias(
                "is_local_max"
            ),
        )
        .drop("middle", "is_middle", "local_max_group", "max_window_score")
        .drop(selectors.ends_with("_shuffled"))
    )
    fdr_indices = find_nearest(
        fdr_table.threshold.values,
        pileup["score"].to_numpy(),
    )
    false_discovery_rates = fdr_table.FDR.values[fdr_indices]
    pileup = (
        pileup.with_columns(FDR=false_discovery_rates)
        .with_columns(
            pl.when(pl.col("score") >= 0.0)
            .then(pl.col("FDR"))
            .otherwise(1.0)
            .alias("FDR"),
        )
        .with_columns(
            (-10 * np.log10(pl.col("FDR")))
            .alias("log_FDR")
            .replace(float("inf"), 100),
        )
    )
    logging.info("Done calculating max window score:\n%s", pileup)
    pileup.write_csv(outfile, separator="\t")


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", type=Path)
    parser.add_argument("outfile", type=Path)
    parser.add_argument("--fdr-table", type=Path)
    parser.add_argument("--nrows", type=int)
    parser.add_argument("--max-cov", type=int)
    parser.add_argument("--min-cov", type=int)
    parser.add_argument("--max-fdr", type=float, default=0.05)
    parser.add_argument("-v", "--verbose", type=int, default=0)
    return parser.parse_args()


def configure_logging(verbose):
    log_format = "[%(levelname)s][Time elapsed (ms) %(relativeCreated)d]: %(message)s"
    log_level = 10 * (3 - verbose)
    logging.basicConfig(format=log_format)
    logging.getLogger().setLevel(log_level)


def main():
    arguments = parse_arguments()
    configure_logging(arguments.verbose)
    if arguments.fdr_table is not None:
        fdr_table = read_fdr_table(arguments.fdr_table)
        apply_fdr_table(
            arguments.infile,
            arguments.outfile,
            fdr_table,
            arguments.nrows,
        )
    else:
        make_fdr_table(
            arguments.infile,
            arguments.outfile,
            arguments.nrows,
            min_coverage=arguments.min_cov,
            max_coverage=arguments.max_cov,
            max_fdr=arguments.max_fdr,
        )


if __name__ == "__main__":
    main()
