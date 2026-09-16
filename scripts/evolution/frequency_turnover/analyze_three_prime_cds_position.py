#!/usr/bin/env python3
"""Compare recent-loss and fixed intron positions using CDS distance from 3' ends.

The analysis follows the positional measurements used by Roy and Gilbert (2005):

1. coding-sequence distance from an intron to the 3' end, expressed in codons;
2. the same distance as a percentage of total coding-sequence length.

For both measurements, smaller values indicate positions closer to the 3' end.
Group 1 recent-loss loci are present in 8--10 of 11 samples, whereas fixed loci
are present in all 11 samples. Introner loci present in either Group 2 sample are
excluded by default, matching the existing turnover analysis.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

MPLCONFIG = Path("/scratch1/chris/tmp/matplotlib")
MPLCONFIG.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu


PROJECT_ROOT = Path(__file__).resolve().parents[3]
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

from figure_color_guide import get_color, load_color_guide


GROUP1_SAMPLES = [
    "CCMP1545",
    "RCC114",
    "RCC1614",
    "RCC1698",
    "RCC2482",
    "RCC373",
    "RCC465",
    "RCC629",
    "RCC692",
    "RCC693",
    "RCC833",
]
GROUP2_SAMPLES = ["RCC1749", "RCC3052"]
CLASS_ORDER = ["recent_loss", "fixed"]
CLASS_DISPLAY = {"recent_loss": "Recent loss", "fixed": "Fixed"}


def parse_args() -> argparse.Namespace:
    outdir = PROJECT_ROOT / "analysis" / "three_prime_cds_position_test"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gtf",
        type=Path,
        default=PROJECT_ROOT / "results" / "annotations" / "CCMP1545.gtf",
    )
    parser.add_argument(
        "--introner-matrix",
        type=Path,
        default=PROJECT_ROOT / "results" / "genotyping" / "genotype_matrix.final.tsv",
    )
    parser.add_argument(
        "--non-introner-matrix",
        type=Path,
        default=PROJECT_ROOT
        / "results"
        / "evolution"
        / "non_introner_introns"
        / "intron_genotype_matrix.tsv",
    )
    parser.add_argument("--outdir", type=Path, default=outdir)
    parser.add_argument("--reference-sample", default="CCMP1545")
    parser.add_argument("--introner-flank-size", type=int, default=100)
    parser.add_argument("--accepted-within-statuses", default="consistent,singleton")
    parser.add_argument(
        "--allow-group2-present-introners",
        action="store_true",
        help="Do not require introner loci to be absent from both Group 2 samples.",
    )
    parser.add_argument("--mating-contig", default="CCMP1545#0#scaffold_2")
    parser.add_argument("--mating-start", type=int, default=49808)
    parser.add_argument("--mating-end", type=int, default=1730591)
    parser.add_argument(
        "--color-guide",
        type=Path,
        default=PROJECT_ROOT / "master_figure_color_guide.tsv",
    )
    return parser.parse_args()


def parse_attrs(value: str) -> dict[str, str]:
    return dict(re.findall(r'(\S+) "([^"]+)"', value))


def parse_csv(value: str) -> set[str]:
    return {item.strip() for item in value.split(",") if item.strip()}


def classify_present_count(present_count: int) -> str | None:
    if present_count == len(GROUP1_SAMPLES):
        return "fixed"
    if present_count in {8, 9, 10}:
        return "recent_loss"
    return None


def in_mating_type_region(
    contig: str,
    position: int,
    mating_contig: str,
    mating_start: int,
    mating_end: int,
) -> bool:
    return contig == mating_contig and mating_start <= position <= mating_end


def cds_bases_toward_3prime(
    intron_start: int,
    intron_end: int,
    strand: str,
    cds_segments: list[tuple[int, int]],
) -> int:
    """Return coding bases after the intron in transcript orientation."""
    if strand == "-":
        return sum(
            end - start + 1
            for start, end in cds_segments
            if end < intron_start
        )
    return sum(
        end - start + 1
        for start, end in cds_segments
        if start > intron_end
    )


def extract_gtf_introns(gtf_path: Path) -> pd.DataFrame:
    transcript_exons: dict[str, list[tuple[int, int]]] = defaultdict(list)
    transcript_cds: dict[str, list[tuple[int, int]]] = defaultdict(list)
    transcript_meta: dict[str, dict[str, str]] = {}

    with open(gtf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in {"exon", "CDS"}:
                continue
            attrs = parse_attrs(fields[8])
            transcript_id = attrs.get("transcript_id")
            if not transcript_id:
                continue
            transcript_meta.setdefault(
                transcript_id,
                {
                    "transcript_id": transcript_id,
                    "gene_id": attrs.get("gene_id", ""),
                    "contig": fields[0],
                    "strand": fields[6],
                },
            )
            interval = (int(fields[3]), int(fields[4]))
            if fields[2] == "exon":
                transcript_exons[transcript_id].append(interval)
            else:
                transcript_cds[transcript_id].append(interval)

    records: list[dict[str, object]] = []
    for transcript_id, exons in transcript_exons.items():
        if len(exons) < 2:
            continue
        cds_segments = sorted(transcript_cds.get(transcript_id, []))
        if not cds_segments:
            continue
        cds_length_nt = sum(end - start + 1 for start, end in cds_segments)
        if cds_length_nt <= 0:
            continue

        genomic_introns = []
        exons = sorted(exons)
        for left, right in zip(exons[:-1], exons[1:], strict=True):
            intron_start = left[1] + 1
            intron_end = right[0] - 1
            if intron_end >= intron_start:
                genomic_introns.append((intron_start, intron_end))
        if not genomic_introns:
            continue

        meta = transcript_meta[transcript_id]
        transcript_introns = (
            list(reversed(genomic_introns))
            if meta["strand"] == "-"
            else genomic_introns
        )
        for intron_index, (intron_start, intron_end) in enumerate(transcript_introns):
            distance_nt = cds_bases_toward_3prime(
                intron_start,
                intron_end,
                meta["strand"],
                cds_segments,
            )
            records.append(
                {
                    **meta,
                    "intron_index": intron_index,
                    "n_introns": len(transcript_introns),
                    "intron_start": intron_start,
                    "intron_end": intron_end,
                    "intron_length_bp": intron_end - intron_start + 1,
                    "cds_length_nt": cds_length_nt,
                    "cds_length_codons": cds_length_nt / 3.0,
                    "cds_nt_from_3prime": distance_nt,
                    "codons_from_3prime": distance_nt / 3.0,
                    "percent_cds_from_3prime": 100.0 * distance_nt / cds_length_nt,
                }
            )

    introns = pd.DataFrame(records)
    if introns.empty:
        raise ValueError(f"No coding transcripts with introns found in {gtf_path}")
    return introns


def best_coordinate_match(
    candidates: pd.DataFrame | None,
    query_start: int,
    query_end: int,
    query_center: int,
) -> dict[str, object] | None:
    if candidates is None or candidates.empty:
        return None
    overlaps = candidates[
        (candidates["intron_start"] <= query_end)
        & (candidates["intron_end"] >= query_start)
    ].copy()
    if overlaps.empty:
        return None
    overlaps["overlap_bp"] = (
        np.minimum(overlaps["intron_end"], query_end)
        - np.maximum(overlaps["intron_start"], query_start)
        + 1
    ).clip(lower=0)
    overlaps["center_contained"] = (
        (overlaps["intron_start"] <= query_center)
        & (overlaps["intron_end"] >= query_center)
    )
    overlaps = overlaps.sort_values(
        ["center_contained", "overlap_bp", "intron_length_bp"],
        ascending=[False, False, True],
    )
    return overlaps.iloc[0].to_dict()


def mapped_record(
    base: dict[str, object],
    match: dict[str, object] | None,
    query_start: int,
    query_end: int,
) -> dict[str, object]:
    record = dict(base)
    record["query_start"] = query_start
    record["query_end"] = query_end
    if match is None:
        record["mapping_status"] = "unmapped"
        for field in [
            "transcript_id",
            "gene_id",
            "strand",
            "intron_index",
            "n_introns",
            "intron_start",
            "intron_end",
            "intron_length_bp",
            "cds_length_nt",
            "cds_length_codons",
            "cds_nt_from_3prime",
            "codons_from_3prime",
            "percent_cds_from_3prime",
            "overlap_bp",
        ]:
            record[field] = np.nan
        return record

    record["mapping_status"] = "mapped"
    for field in [
        "transcript_id",
        "gene_id",
        "strand",
        "intron_index",
        "n_introns",
        "intron_start",
        "intron_end",
        "intron_length_bp",
        "cds_length_nt",
        "cds_length_codons",
        "cds_nt_from_3prime",
        "codons_from_3prime",
        "percent_cds_from_3prime",
    ]:
        record[field] = match[field]
    record["overlap_bp"] = max(
        0,
        min(int(match["intron_end"]), query_end)
        - max(int(match["intron_start"]), query_start)
        + 1,
    )
    return record


def load_introner_loci(
    args: argparse.Namespace,
    gtf_by_contig: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    usecols = [
        "ortholog_id",
        "sample",
        "start",
        "end",
        "contig",
        "presence",
        "within_group_status",
    ]
    matrix = pd.read_csv(args.introner_matrix, sep="\t", usecols=usecols)
    accepted_statuses = parse_csv(args.accepted_within_statuses)
    records: list[dict[str, object]] = []

    for ortholog_id, group in matrix.groupby("ortholog_id", sort=True):
        g1 = group[group["sample"].isin(GROUP1_SAMPLES)]
        by_sample = g1.drop_duplicates("sample").set_index("sample")
        if not set(GROUP1_SAMPLES).issubset(by_sample.index):
            continue
        calls = by_sample.loc[GROUP1_SAMPLES, "presence"]
        if not calls.isin([1, 2]).all():
            continue
        present_count = int((calls == 1).sum())
        class_label = classify_present_count(present_count)
        if class_label is None:
            continue

        statuses = {
            str(value)
            for value in group["within_group_status"].dropna().unique()
            if str(value)
        }
        if statuses and not statuses.intersection(accepted_statuses):
            continue
        group2_present_count = int(
            (
                group[group["sample"].isin(GROUP2_SAMPLES)]["presence"]
                == 1
            ).sum()
        )
        if group2_present_count and not args.allow_group2_present_introners:
            continue

        ref = group[group["sample"] == args.reference_sample]
        if ref.empty:
            continue
        ref = ref.iloc[0]
        start, end = int(ref["start"]), int(ref["end"])
        if end < start:
            start, end = end, start
        contig = str(ref["contig"])
        center = (start + end) // 2
        if in_mating_type_region(
            contig,
            center,
            args.mating_contig,
            args.mating_start,
            args.mating_end,
        ):
            continue

        body_start = start + args.introner_flank_size
        body_end = end - args.introner_flank_size
        if body_end < body_start:
            body_start = body_end = center
        match = best_coordinate_match(
            gtf_by_contig.get(contig),
            body_start,
            body_end,
            center,
        )
        base = {
            "locus_type": "introner",
            "locus_id": str(ortholog_id),
            "class_label": class_label,
            "class_display": CLASS_DISPLAY[class_label],
            "contig": contig,
            "locus_start": start,
            "locus_end": end,
            "locus_position": center,
            "group1_present_count": present_count,
            "group1_absent_count": len(GROUP1_SAMPLES) - present_count,
            "group2_present_count": group2_present_count,
        }
        records.append(mapped_record(base, match, body_start, body_end))
    return pd.DataFrame(records)


def load_regular_intron_loci(
    args: argparse.Namespace,
    gtf_introns: pd.DataFrame,
) -> pd.DataFrame:
    matrix = pd.read_csv(args.non_introner_matrix, sep="\t")
    gtf_by_key = {
        (str(row.transcript_id), int(row.intron_index)): row._asdict()
        for row in gtf_introns.itertuples(index=False)
    }
    records: list[dict[str, object]] = []

    for row in matrix.itertuples(index=False):
        values = row._asdict()
        calls = [values[sample] for sample in GROUP1_SAMPLES]
        if not all(value in {1, 2} for value in calls):
            continue
        present_count = int(sum(value == 1 for value in calls))
        class_label = classify_present_count(present_count)
        if class_label is None:
            continue
        start, end = int(values["ref_start"]), int(values["ref_end"])
        if end < start:
            start, end = end, start
        contig = str(values["contig"])
        center = (start + end) // 2
        if in_mating_type_region(
            contig,
            center,
            args.mating_contig,
            args.mating_start,
            args.mating_end,
        ):
            continue

        transcript_id = str(values["gene_id"])
        intron_index = int(values["intron_index"])
        match = gtf_by_key.get((transcript_id, intron_index))
        base = {
            "locus_type": "regular_intron",
            "locus_id": f"{transcript_id}:{intron_index}",
            "class_label": class_label,
            "class_display": CLASS_DISPLAY[class_label],
            "contig": contig,
            "locus_start": start,
            "locus_end": end,
            "locus_position": center,
            "group1_present_count": present_count,
            "group1_absent_count": len(GROUP1_SAMPLES) - present_count,
            "group2_present_count": int(
                sum(values.get(sample, np.nan) == 1 for sample in GROUP2_SAMPLES)
            ),
        }
        records.append(mapped_record(base, match, start, end))
    return pd.DataFrame(records)


def summarize_loci(loci: pd.DataFrame) -> pd.DataFrame:
    records = []
    for (locus_type, class_label), group in loci.groupby(
        ["locus_type", "class_label"], sort=False
    ):
        mapped = group[group["mapping_status"] == "mapped"]
        records.append(
            {
                "locus_type": locus_type,
                "class_label": class_label,
                "class_display": CLASS_DISPLAY[class_label],
                "n_loci": int(group["locus_id"].nunique()),
                "n_mapped": int(mapped["locus_id"].nunique()),
                "mapping_fraction": len(mapped) / len(group) if len(group) else np.nan,
                "mean_codons_from_3prime": mapped["codons_from_3prime"].mean(),
                "median_codons_from_3prime": mapped["codons_from_3prime"].median(),
                "mean_percent_cds_from_3prime": mapped[
                    "percent_cds_from_3prime"
                ].mean(),
                "median_percent_cds_from_3prime": mapped[
                    "percent_cds_from_3prime"
                ].median(),
            }
        )
    result = pd.DataFrame(records)
    result["class_order"] = result["class_label"].map(
        {label: index for index, label in enumerate(CLASS_ORDER)}
    )
    return result.sort_values(["locus_type", "class_order"]).drop(
        columns="class_order"
    )


def run_tests(loci: pd.DataFrame) -> pd.DataFrame:
    records = []
    metrics = [
        ("codons_from_3prime", "Codons from 3' end"),
        ("percent_cds_from_3prime", "Percent CDS from 3' end"),
    ]
    for locus_type in ["introner", "regular_intron"]:
        subset = loci[
            (loci["locus_type"] == locus_type)
            & (loci["mapping_status"] == "mapped")
        ]
        for metric, metric_display in metrics:
            recent = pd.to_numeric(
                subset.loc[subset["class_label"] == "recent_loss", metric],
                errors="coerce",
            ).dropna()
            fixed = pd.to_numeric(
                subset.loc[subset["class_label"] == "fixed", metric],
                errors="coerce",
            ).dropna()
            if recent.empty or fixed.empty:
                u_stat = p_two_sided = p_loss_closer = np.nan
            else:
                two_sided = mannwhitneyu(
                    recent,
                    fixed,
                    alternative="two-sided",
                    method="asymptotic",
                )
                closer = mannwhitneyu(
                    recent,
                    fixed,
                    alternative="less",
                    method="asymptotic",
                )
                u_stat = float(two_sided.statistic)
                p_two_sided = float(two_sided.pvalue)
                p_loss_closer = float(closer.pvalue)
            records.append(
                {
                    "locus_type": locus_type,
                    "comparison": "recent_loss_vs_fixed",
                    "metric": metric,
                    "metric_display": metric_display,
                    "recent_loss_n": len(recent),
                    "fixed_n": len(fixed),
                    "recent_loss_median": recent.median() if len(recent) else np.nan,
                    "fixed_median": fixed.median() if len(fixed) else np.nan,
                    "median_difference_recent_loss_minus_fixed": (
                        recent.median() - fixed.median()
                        if len(recent) and len(fixed)
                        else np.nan
                    ),
                    "mannwhitney_u": u_stat,
                    "mannwhitney_p_two_sided": p_two_sided,
                    "mannwhitney_p_recent_loss_closer_3prime": p_loss_closer,
                }
            )
    return pd.DataFrame(records)


def validate(loci: pd.DataFrame) -> pd.DataFrame:
    mapped = loci[loci["mapping_status"] == "mapped"]
    checks = [
        (
            "all_percent_positions_within_0_100",
            bool(mapped["percent_cds_from_3prime"].between(0, 100).all()),
        ),
        (
            "all_codon_distances_nonnegative",
            bool((mapped["codons_from_3prime"] >= 0).all()),
        ),
        (
            "all_recent_loss_counts_are_8_to_10",
            bool(
                loci.loc[loci["class_label"] == "recent_loss", "group1_present_count"]
                .isin([8, 9, 10])
                .all()
            ),
        ),
        (
            "all_fixed_counts_are_11",
            bool(
                (
                    loci.loc[loci["class_label"] == "fixed", "group1_present_count"]
                    == 11
                ).all()
            ),
        ),
        (
            "all_requested_sets_have_mapped_loci",
            bool(
                (
                    mapped.groupby(["locus_type", "class_label"])["locus_id"]
                    .nunique()
                    .reindex(
                        pd.MultiIndex.from_product(
                            [["introner", "regular_intron"], CLASS_ORDER]
                        ),
                        fill_value=0,
                    )
                    > 0
                ).all()
            ),
        ),
    ]
    return pd.DataFrame(checks, columns=["check", "pass"])


def plot_distributions(
    loci: pd.DataFrame,
    tests: pd.DataFrame,
    color_guide: Path,
    output: Path,
) -> None:
    plot_data = loci[loci["mapping_status"] == "mapped"].copy()
    colors = load_color_guide(color_guide)
    palette = {
        "Recent loss": get_color("Population 1 Polymorphic", colors),
        "Fixed": get_color("Population 1 Fixed", colors),
    }
    sns.set_style("whitegrid")
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    metrics = [
        ("codons_from_3prime", "Codons from 3′ end"),
        ("percent_cds_from_3prime", "Coding sequence from 3′ end (%)"),
    ]
    locus_types = [
        ("introner", "Introners"),
        ("regular_intron", "Regular introns"),
    ]

    for row_index, (locus_type, row_title) in enumerate(locus_types):
        subset = plot_data[plot_data["locus_type"] == locus_type]
        for col_index, (metric, ylabel) in enumerate(metrics):
            ax = axes[row_index, col_index]
            sns.boxplot(
                data=subset,
                x="class_display",
                y=metric,
                order=["Recent loss", "Fixed"],
                hue="class_display",
                hue_order=["Recent loss", "Fixed"],
                palette=palette,
                showfliers=False,
                width=0.55,
                linewidth=1.2,
                legend=False,
                ax=ax,
            )
            test = tests[
                (tests["locus_type"] == locus_type)
                & (tests["metric"] == metric)
            ].iloc[0]
            ax.set_title(
                f"{row_title}: MWU two-sided p={test['mannwhitney_p_two_sided']:.3g}",
                fontsize=11,
            )
            ax.set_xlabel("")
            ax.set_ylabel(ylabel)
            ax.grid(axis="x", visible=False)
            if metric == "percent_cds_from_3prime":
                ax.set_ylim(0, 100)

    fig.suptitle(
        "Group 1 recent-loss versus fixed positions along coding sequence",
        fontsize=14,
        y=1.01,
    )
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    gtf_introns = extract_gtf_introns(args.gtf)
    gtf_by_contig = {
        contig: group.reset_index(drop=True)
        for contig, group in gtf_introns.groupby("contig")
    }
    introners = load_introner_loci(args, gtf_by_contig)
    regular_introns = load_regular_intron_loci(args, gtf_introns)
    loci = pd.concat([introners, regular_introns], ignore_index=True)
    loci["class_order"] = loci["class_label"].map(
        {label: index for index, label in enumerate(CLASS_ORDER)}
    )
    loci = loci.sort_values(
        ["locus_type", "class_order", "contig", "locus_position", "locus_id"]
    ).drop(columns="class_order")

    summary = summarize_loci(loci)
    tests = run_tests(loci)
    validation = validate(loci)
    if not validation["pass"].all():
        failures = validation.loc[~validation["pass"], "check"].tolist()
        raise ValueError(f"Validation failed: {failures}")

    gtf_introns.to_csv(args.outdir / "ccmp1545_cds_intron_positions.tsv", sep="\t", index=False)
    loci.to_csv(args.outdir / "three_prime_cds_locus_positions.tsv", sep="\t", index=False)
    summary.to_csv(args.outdir / "three_prime_cds_position_summary.tsv", sep="\t", index=False)
    tests.to_csv(args.outdir / "three_prime_cds_mannwhitney.tsv", sep="\t", index=False)
    validation.to_csv(args.outdir / "validation_summary.tsv", sep="\t", index=False)
    plot_distributions(
        loci,
        tests,
        args.color_guide,
        args.outdir / "three_prime_cds_position_distributions.png",
    )

    print("Summary by class:")
    print(summary.to_string(index=False))
    print("\nMann-Whitney U tests:")
    print(tests.to_string(index=False))
    print("\nValidation:")
    print(validation.to_string(index=False))


if __name__ == "__main__":
    main()
