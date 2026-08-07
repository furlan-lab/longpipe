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


def parse_count(value):
    number = parse_number(value)
    if number is None or number < 0 or not number.is_integer():
        return None
    return int(number)


def log_combination(total, selected):
    if selected < 0 or selected > total:
        return float("-inf")
    return (
        math.lgamma(total + 1)
        - math.lgamma(selected + 1)
        - math.lgamma(total - selected + 1)
    )


def fisher_exact_two_sided(
    accessible_h1,
    inaccessible_h1,
    accessible_h2,
    inaccessible_h2,
):
    first_row = accessible_h1 + inaccessible_h1
    second_row = accessible_h2 + inaccessible_h2
    accessible_total = accessible_h1 + accessible_h2
    total = first_row + second_row
    minimum = max(0, first_row - (total - accessible_total))
    maximum = min(first_row, accessible_total)

    def log_probability(first_accessible):
        return (
            log_combination(accessible_total, first_accessible)
            + log_combination(
                total - accessible_total,
                first_row - first_accessible,
            )
            - log_combination(total, first_row)
        )

    probabilities = [
        log_probability(first_accessible)
        for first_accessible in range(minimum, maximum + 1)
    ]
    observed = log_probability(accessible_h1)
    maximum_log_probability = max(probabilities)
    weights = [
        math.exp(probability - maximum_log_probability)
        for probability in probabilities
    ]
    selected_weight = math.fsum(
        weight
        for probability, weight in zip(probabilities, weights)
        if probability <= observed + 1e-7
    )
    return min(1.0, selected_weight / math.fsum(weights))


def odds_ratio(accessible_h1, inaccessible_h1, accessible_h2, inaccessible_h2):
    numerator = accessible_h1 * inaccessible_h2
    denominator = inaccessible_h1 * accessible_h2
    if denominator == 0:
        return float("nan") if numerator == 0 else float("inf")
    return numerator / denominator


def benjamini_hochberg(pvalues):
    adjusted = [1.0] * len(pvalues)
    ranked = sorted(enumerate(pvalues), key=lambda item: item[1])
    running_minimum = 1.0
    test_count = len(ranked)
    for rank, (index, pvalue) in reversed(list(enumerate(ranked, start=1))):
        running_minimum = min(
            running_minimum,
            pvalue * test_count / rank,
        )
        adjusted[index] = min(1.0, running_minimum)
    return adjusted


def is_true(value):
    return str(value).strip().lower() in {"true", "t", "yes", "y", "1"}


def format_number(value):
    if math.isnan(value):
        return "nan"
    if math.isinf(value):
        return "inf" if value > 0 else "-inf"
    return f"{value:.10g}"


def analyze_pileup(
    pileup_path,
    output_path,
    summary_path,
    minimum_coverage,
    maximum_fire_fdr,
    qvalue_threshold,
):
    required = {
        "#chrom",
        "start",
        "end",
        "FDR",
        "is_local_max",
        "coverage_H1",
        "fire_coverage_H1",
        "coverage_H2",
        "fire_coverage_H2",
    }
    records = []
    candidate_count = 0
    with open_text(pileup_path) as input_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError(
                "FIRE pileup is missing haplotype-analysis columns: "
                + ", ".join(missing)
            )
        for row in reader:
            fire_fdr = parse_number(row["FDR"])
            if (
                not is_true(row["is_local_max"])
                or fire_fdr is None
                or fire_fdr > maximum_fire_fdr
            ):
                continue
            candidate_count += 1
            coverage_h1 = parse_count(row["coverage_H1"])
            fire_h1 = parse_count(row["fire_coverage_H1"])
            coverage_h2 = parse_count(row["coverage_H2"])
            fire_h2 = parse_count(row["fire_coverage_H2"])
            if None in {coverage_h1, fire_h1, coverage_h2, fire_h2}:
                continue
            if coverage_h1 < minimum_coverage or coverage_h2 < minimum_coverage:
                continue
            inaccessible_h1 = coverage_h1 - fire_h1
            inaccessible_h2 = coverage_h2 - fire_h2
            if inaccessible_h1 < 0 or inaccessible_h2 < 0:
                continue
            pvalue = fisher_exact_two_sided(
                fire_h1,
                inaccessible_h1,
                fire_h2,
                inaccessible_h2,
            )
            percent_h1 = 100.0 * fire_h1 / coverage_h1
            percent_h2 = 100.0 * fire_h2 / coverage_h2
            records.append(
                {
                    "chromosome": row["#chrom"],
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "coverage_h1": coverage_h1,
                    "fire_h1": fire_h1,
                    "percent_h1": percent_h1,
                    "coverage_h2": coverage_h2,
                    "fire_h2": fire_h2,
                    "percent_h2": percent_h2,
                    "difference": percent_h1 - percent_h2,
                    "odds_ratio": odds_ratio(
                        fire_h1,
                        inaccessible_h1,
                        fire_h2,
                        inaccessible_h2,
                    ),
                    "pvalue": pvalue,
                    "fire_fdr": fire_fdr,
                }
            )

    qvalues = benjamini_hochberg(
        [record["pvalue"] for record in records]
    )
    header = (
        "#chrom",
        "start",
        "end",
        "coverage_H1",
        "fire_coverage_H1",
        "percent_accessible_H1",
        "coverage_H2",
        "fire_coverage_H2",
        "percent_accessible_H2",
        "H1_minus_H2_percent_accessible",
        "odds_ratio",
        "p_value",
        "q_value",
        "significant",
        "FIRE_FDR",
    )
    significant_count = 0
    with output_path.open("w", encoding="utf-8", newline="") as output_handle:
        writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for record, qvalue in zip(records, qvalues):
            significant = qvalue <= qvalue_threshold
            significant_count += int(significant)
            writer.writerow(
                (
                    record["chromosome"],
                    record["start"],
                    record["end"],
                    record["coverage_h1"],
                    record["fire_h1"],
                    format_number(record["percent_h1"]),
                    record["coverage_h2"],
                    record["fire_h2"],
                    format_number(record["percent_h2"]),
                    format_number(record["difference"]),
                    format_number(record["odds_ratio"]),
                    format_number(record["pvalue"]),
                    format_number(qvalue),
                    "yes" if significant else "no",
                    format_number(record["fire_fdr"]),
                )
            )

    with summary_path.open("w", encoding="utf-8", newline="") as summary_handle:
        writer = csv.writer(summary_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "value"))
        writer.writerow(("candidate_local_maxima", candidate_count))
        writer.writerow(("tested_loci", len(records)))
        writer.writerow(("significant_loci", significant_count))
        writer.writerow(("minimum_haplotype_coverage", minimum_coverage))
        writer.writerow(("maximum_FIRE_FDR", maximum_fire_fdr))
        writer.writerow(("significant_q_value", qvalue_threshold))
        writer.writerow(("test", "two-sided Fisher exact with BH correction"))


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("pileup", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("summary", type=Path)
    parser.add_argument("--minimum-coverage", type=int, required=True)
    parser.add_argument("--maximum-fire-fdr", type=float, required=True)
    parser.add_argument("--qvalue-threshold", type=float, required=True)
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    analyze_pileup(
        arguments.pileup,
        arguments.output,
        arguments.summary,
        arguments.minimum_coverage,
        arguments.maximum_fire_fdr,
        arguments.qvalue_threshold,
    )


if __name__ == "__main__":
    main()
