import argparse
import math
import sys
from pathlib import Path

import polars as pl


def weighted_median(dataframe, value_column, weight_column):
    grouped = (
        dataframe.groupby(value_column)[weight_column]
        .sum()
        .reset_index()
        .sort_values(value_column)
    )
    print(grouped, file=sys.stderr)
    grouped["cumulative_weight"] = grouped[weight_column].cumsum()
    grouped["cutoff"] = grouped[weight_column].sum() / 2.0
    print(grouped, file=sys.stderr)
    values = grouped.loc[
        grouped["cumulative_weight"] >= grouped["cutoff"], value_column
    ]
    if values.empty:
        raise ValueError("FIRE coverage input contains no usable intervals")
    return values.iloc[0]


def read_coverage(path, chromosomes):
    return (
        pl.read_csv(
            path,
            separator="\t",
            has_header=False,
            new_columns=["chromosome", "start", "end", "coverage"],
            low_memory=True,
            schema_overrides=[pl.Utf8, pl.Int64, pl.Int64, pl.Int64],
        )
        .lazy()
        .filter(pl.col("coverage") > 0)
        .filter(pl.col("chromosome").is_in(chromosomes))
        .drop("chromosome")
        .with_columns((pl.col("end") - pl.col("start")).alias("weight"))
        .drop(["start", "end"])
        .collect()
        .to_pandas()
    )


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("bedgraph", type=Path)
    parser.add_argument("median_output", type=Path)
    parser.add_argument("minimum_output", type=Path)
    parser.add_argument("maximum_output", type=Path)
    parser.add_argument("--chromosome", action="append", required=True)
    parser.add_argument("--coverage-within-n-sd", type=float, required=True)
    parser.add_argument("--min-coverage", type=int, required=True)
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    dataframe = read_coverage(arguments.bedgraph, arguments.chromosome)
    print(dataframe, file=sys.stderr)
    median_coverage = weighted_median(dataframe, "coverage", "weight")
    standard_deviation = math.sqrt(median_coverage)
    minimum_coverage = max(
        median_coverage
        - arguments.coverage_within_n_sd * standard_deviation,
        arguments.min_coverage,
    )
    maximum_coverage = (
        median_coverage
        + arguments.coverage_within_n_sd * standard_deviation
    )
    mean_coverage = (
        (dataframe["coverage"] * dataframe["weight"]).sum()
        / dataframe["weight"].sum()
    )
    print(f"\nmean coverage: {mean_coverage}", file=sys.stderr)
    print(f"median coverage: {median_coverage}\n", file=sys.stderr)
    if median_coverage < 5:
        raise ValueError(
            f"Median coverage is {median_coverage}; FIRE requires at least 5x."
        )
    arguments.median_output.write_text(
        f"{round(median_coverage)}\n", encoding="utf-8"
    )
    arguments.minimum_output.write_text(
        f"{round(minimum_coverage)}\n", encoding="utf-8"
    )
    arguments.maximum_output.write_text(
        f"{round(maximum_coverage)}\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
