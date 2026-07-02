#!/usr/bin/env python3
"""Test introner enrichment/depletion in SweepFinder2 sweep-like regions."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

Path("/scratch1/chris/tmp/matplotlib").mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", "/scratch1/chris/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clr-tsv", required=True, help="Merged SweepFinder2 CLR table.")
    parser.add_argument("--sweep-regions", required=True, help="Called sweep regions BED.")
    parser.add_argument("--contigs-tsv", required=True, help="Contig table with lengths.")
    parser.add_argument("--introner-matrix", required=True, help="genotype_matrix.final.tsv")
    parser.add_argument("--non-introner-matrix", required=True, help="intron_genotype_matrix.tsv")
    parser.add_argument("--group1-samples", required=True, help="Comma-separated Group 1 sample names.")
    parser.add_argument("--group2-samples", required=True, help="Comma-separated Group 2 sample names.")
    parser.add_argument(
        "--accepted-within-statuses",
        default="consistent,singleton",
        help="Comma-separated within_group_status values to retain when status is present.",
    )
    parser.add_argument("--reference-sample", default="CCMP1545")
    parser.add_argument("--flank-size", type=int, default=100, help="Flanking bp to trim for introner midpoint.")
    parser.add_argument("--mating-contig", required=True)
    parser.add_argument("--mating-start", type=int, required=True)
    parser.add_argument("--mating-end", type=int, required=True)
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--region-overlap-tsv", required=True)
    parser.add_argument("--perlocus-tsv", required=True)
    parser.add_argument("--summary-tsv", required=True)
    parser.add_argument("--plot-enrichment-png", required=True)
    parser.add_argument("--plot-enrichment-pdf", required=True)
    parser.add_argument("--plot-clr-png", required=True)
    parser.add_argument("--plot-clr-pdf", required=True)
    return parser.parse_args()


def parse_csv(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def overlaps_interval(contig: str, pos: int, query_contig: str, start: int, end: int) -> bool:
    return contig == query_contig and start <= pos <= end


def load_callable_segments(
    contigs: pd.DataFrame,
    mating_contig: str,
    mating_start: int,
    mating_end: int,
) -> pd.DataFrame:
    rows = []
    for row in contigs.itertuples(index=False):
        contig = str(row.contig)
        length = int(row.length)
        if contig == mating_contig:
            if mating_start > 1:
                rows.append({"contig": contig, "start": 1, "end": mating_start - 1})
            if mating_end < length:
                rows.append({"contig": contig, "start": mating_end + 1, "end": length})
        else:
            rows.append({"contig": contig, "start": 1, "end": length})
    return pd.DataFrame(rows)


def fraction_overlapping(loci: pd.DataFrame, regions: pd.DataFrame) -> float:
    if loci.empty:
        return np.nan
    if regions.empty:
        return 0.0
    flags = []
    for contig, group in loci.groupby("contig", sort=False):
        positions = group["position"].to_numpy(dtype=np.int64)
        contig_regions = regions[regions["contig"] == contig]
        if contig_regions.empty:
            flags.append(np.zeros(len(positions), dtype=bool))
            continue
        overlap = np.zeros(len(positions), dtype=bool)
        for row in contig_regions.itertuples(index=False):
            overlap |= (positions >= int(row.start)) & (positions <= int(row.end))
        flags.append(overlap)
    return float(np.concatenate(flags).mean())


def sample_positions(segments: pd.DataFrame, n: int, rng: np.random.Generator) -> pd.DataFrame:
    if segments.empty or n <= 0:
        return pd.DataFrame(columns=["contig", "position"])
    lengths = (segments["end"] - segments["start"] + 1).to_numpy(dtype=np.int64)
    total = int(lengths.sum())
    if total <= 0:
        return pd.DataFrame(columns=["contig", "position"])

    draws = rng.integers(0, total, size=n)
    segment_ends = np.cumsum(lengths)
    segment_starts = segment_ends - lengths
    seg_idx = np.searchsorted(segment_ends, draws, side="right")
    offsets = draws - segment_starts[seg_idx]
    contigs = segments["contig"].to_numpy()[seg_idx]
    positions = segments["start"].to_numpy(dtype=np.int64)[seg_idx] + offsets
    return pd.DataFrame({"contig": contigs, "position": positions})


def load_stratified_introners(
    matrix_path: str,
    group1_samples: list[str],
    group2_samples: list[str],
    reference_sample: str,
    accepted_within_statuses: set[str],
    flank_size: int,
    mating_contig: str,
    mating_start: int,
    mating_end: int,
) -> dict[str, pd.DataFrame]:
    """Return fixed-G1 and polymorphic-G1 introner loci anchored to CCMP1545."""
    usecols = [
        "ortholog_id",
        "sample",
        "contig",
        "start",
        "end",
        "presence",
        "within_group_status",
    ]
    df = pd.read_csv(matrix_path, sep="\t", usecols=usecols)
    if "within_group_status" not in df.columns:
        df["within_group_status"] = ""

    fixed_rows = []
    polymorphic_rows = []

    for ortholog_id, group in df.groupby("ortholog_id", sort=True):
        g1 = group[group["sample"].isin(group1_samples)]
        g1_by_sample = g1.drop_duplicates("sample").set_index("sample")
        if not set(group1_samples).issubset(g1_by_sample.index):
            continue

        g1_presence = g1_by_sample.loc[group1_samples, "presence"]
        if not g1_presence.isin([1, 2]).all():
            continue

        present_count = int((g1_presence == 1).sum())
        absent_count = int((g1_presence == 2).sum())
        if present_count == 0:
            continue

        g2 = group[group["sample"].isin(group2_samples)]
        if int((g2["presence"] == 1).sum()) > 0:
            continue

        within_values = {
            str(value)
            for value in group["within_group_status"].dropna().unique()
            if str(value)
        }
        if within_values and not (within_values & accepted_within_statuses):
            continue

        if present_count == len(group1_samples):
            locus_class = "introner_fixed_g1"
        elif absent_count > 0:
            locus_class = "introner_polymorphic_g1"
        else:
            continue

        ref_rows = group[group["sample"] == reference_sample]
        if ref_rows.empty:
            continue
        ref = ref_rows.iloc[0]
        try:
            start = int(ref["start"])
            end = int(ref["end"])
        except (TypeError, ValueError):
            continue
        contig = str(ref["contig"])
        if start < 0 or end < 0 or not contig:
            continue
        if end < start:
            start, end = end, start

        body_start = start + flank_size
        body_end = end - flank_size
        if body_end <= body_start:
            continue
        position = (body_start + body_end) // 2
        if overlaps_interval(contig, position, mating_contig, mating_start, mating_end):
            continue

        row = {
            "ortholog_id": ortholog_id,
            "contig": contig,
            "start": start,
            "end": end,
            "position": position,
            "group1_present_count": present_count,
            "group1_absent_count": absent_count,
            "reference_presence": int(ref["presence"]) if pd.notna(ref["presence"]) else np.nan,
        }
        if locus_class == "introner_fixed_g1":
            fixed_rows.append(row)
        else:
            polymorphic_rows.append(row)

    return {
        "introner_fixed_g1": pd.DataFrame(fixed_rows).drop_duplicates("ortholog_id"),
        "introner_polymorphic_g1": pd.DataFrame(polymorphic_rows).drop_duplicates("ortholog_id"),
    }


def load_non_introner_introns(
    matrix_path: str,
    mating_contig: str,
    mating_start: int,
    mating_end: int,
) -> pd.DataFrame:
    df = pd.read_csv(matrix_path, sep="\t")
    if "CCMP1545" not in df.columns:
        raise ValueError("Non-introner intron matrix missing CCMP1545 column")
    present = df[df["CCMP1545"] == 1].copy()
    present["contig"] = present["contig"].astype(str)
    present["position"] = ((present["ref_start"] + present["ref_end"]) // 2).astype(int)
    mask = [
        not overlaps_interval(contig, pos, mating_contig, mating_start, mating_end)
        for contig, pos in zip(present["contig"], present["position"], strict=True)
    ]
    present = present[mask]
    return present[["gene_id", "contig", "ref_start", "ref_end", "position"]].drop_duplicates()


def load_sweep_regions(path: str) -> pd.DataFrame:
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame(columns=["contig", "start", "end", "max_clr"])
    df = pd.read_csv(path, sep="\t", header=None, names=["contig", "bed_start", "end", "max_clr"])
    df["start"] = df["bed_start"] + 1
    return df[["contig", "start", "end", "max_clr"]]


def assign_nearest_clr(loci: pd.DataFrame, clr_df: pd.DataFrame) -> pd.Series:
    if loci.empty or clr_df.empty:
        return pd.Series(dtype=float, index=loci.index)
    values = pd.Series(np.nan, index=loci.index, dtype=float)
    for contig, locus_group in loci.groupby("contig", sort=False):
        grid = clr_df[clr_df["contig"] == contig].sort_values("position")
        if grid.empty:
            continue
        left = locus_group[["position"]].sort_values("position")
        merged = pd.merge_asof(
            left,
            grid[["position", "clr"]],
            on="position",
            direction="nearest",
        )
        values.loc[left.index] = merged["clr"].to_numpy()
    return values


def permutation_overlap_enrichment(
    loci: pd.DataFrame,
    regions: pd.DataFrame,
    segments: pd.DataFrame,
    n_perm: int,
    rng: np.random.Generator,
) -> tuple[float, np.ndarray]:
    if loci.empty:
        return np.nan, np.array([])
    observed = fraction_overlapping(loci, regions)
    null = np.array(
        [fraction_overlapping(sample_positions(segments, len(loci), rng), regions) for _ in range(n_perm)],
        dtype=float,
    )
    return observed, null


def two_sided_perm_p(observed: float, null: np.ndarray) -> float:
    if null.size == 0 or np.isnan(observed):
        return np.nan
    upper = (np.sum(null >= observed) + 1) / (null.size + 1)
    lower = (np.sum(null <= observed) + 1) / (null.size + 1)
    return float(min(2 * min(upper, lower), 1.0))


def plot_enrichment(region_df: pd.DataFrame, png: str, pdf: str) -> None:
    fig, ax = plt.subplots(figsize=(9, 5))
    plot_order = ["introner_fixed_g1", "introner_polymorphic_g1", "non_introner_intron"]
    label_map = {
        "introner_fixed_g1": "fixed G1 introner",
        "introner_polymorphic_g1": "polymorphic G1 introner",
        "non_introner_intron": "non-introner intron",
    }
    colors = {"introner_fixed_g1": "#1b7837", "introner_polymorphic_g1": "#5aae61", "non_introner_intron": "#762a83"}
    plot_df = region_df[region_df["locus_class"].isin(plot_order)].copy()
    plot_df["plot_order"] = plot_df["locus_class"].map({name: idx for idx, name in enumerate(plot_order)})
    plot_df = plot_df.sort_values("plot_order")

    if plot_df.empty:
        ax.text(0.5, 0.5, "No enrichment data", ha="center", va="center", transform=ax.transAxes)
    else:
        x = np.arange(len(plot_df))
        bar_colors = [colors[row.locus_class] for row in plot_df.itertuples(index=False)]
        ax.bar(x, plot_df["observed_overlap_fraction"], color=bar_colors, alpha=0.85)
        ax.errorbar(
            x,
            plot_df["observed_overlap_fraction"],
            yerr=[
                plot_df["observed_overlap_fraction"] - plot_df["perm_ci_low"],
                plot_df["perm_ci_high"] - plot_df["observed_overlap_fraction"],
            ],
            fmt="none",
            ecolor="black",
            capsize=4,
        )
        ax.set_xticks(x)
        ax.set_xticklabels([label_map[row] for row in plot_df["locus_class"]], rotation=15)
        ax.set_ylabel("Fraction overlapping sweep regions")
        ax.set_title("Sweep-region overlap vs permutation null")
    fig.tight_layout()
    fig.savefig(png, dpi=250, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)


def plot_clr_distributions(
    locus_clr: dict[str, pd.Series],
    grid_clr: pd.Series,
    png: str,
    pdf: str,
) -> None:
    fig, ax = plt.subplots(figsize=(10, 5))
    label_map = {
        "introner_fixed_g1": "fixed G1 introner",
        "introner_polymorphic_g1": "polymorphic G1 introner",
        "non_introner_intron": "non-introner intron",
        "genome_grid": "genome grid (subsampled)",
    }
    plot_order = ["introner_fixed_g1", "introner_polymorphic_g1", "non_introner_intron", "genome_grid"]
    data = []
    labels = []
    for key in plot_order:
        if key == "genome_grid":
            series = grid_clr
        else:
            series = locus_clr.get(key, pd.Series(dtype=float))
        vals = pd.to_numeric(series, errors="coerce").dropna()
        if vals.empty:
            continue
        if key == "genome_grid" and len(vals) > 5000:
            vals = vals.sample(5000, random_state=42)
        data.append(vals)
        labels.append(label_map[key])
    if data:
        ax.boxplot(data, tick_labels=labels, showfliers=False)
        ax.tick_params(axis="x", labelrotation=20)
    ax.set_ylabel("Nearest-grid SweepFinder2 CLR")
    ax.set_title("Per-locus CLR by introner class vs controls")
    fig.tight_layout()
    fig.savefig(png, dpi=250, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    group1_samples = parse_csv(args.group1_samples)
    group2_samples = parse_csv(args.group2_samples)
    accepted_within_statuses = set(parse_csv(args.accepted_within_statuses))

    for output in [
        args.region_overlap_tsv,
        args.perlocus_tsv,
        args.summary_tsv,
        args.plot_enrichment_png,
        args.plot_enrichment_pdf,
        args.plot_clr_png,
        args.plot_clr_pdf,
    ]:
        Path(output).parent.mkdir(parents=True, exist_ok=True)

    clr_df = pd.read_csv(args.clr_tsv, sep="\t")
    contigs = pd.read_csv(args.contigs_tsv, sep="\t")
    regions = load_sweep_regions(args.sweep_regions)
    segments = load_callable_segments(contigs, args.mating_contig, args.mating_start, args.mating_end)

    introner_sets = load_stratified_introners(
        args.introner_matrix,
        group1_samples,
        group2_samples,
        args.reference_sample,
        accepted_within_statuses,
        args.flank_size,
        args.mating_contig,
        args.mating_start,
        args.mating_end,
    )
    introns = load_non_introner_introns(
        args.non_introner_matrix,
        args.mating_contig,
        args.mating_start,
        args.mating_end,
    )

    region_rows = []
    locus_classes = [
        ("introner_fixed_g1", introner_sets["introner_fixed_g1"]),
        ("introner_polymorphic_g1", introner_sets["introner_polymorphic_g1"]),
        ("non_introner_intron", introns),
    ]
    for label, loci in locus_classes:
        observed, null = permutation_overlap_enrichment(loci, regions, segments, args.permutations, rng)
        p_value = two_sided_perm_p(observed, null)
        ci_low = float(np.percentile(null, 2.5)) if null.size else np.nan
        ci_high = float(np.percentile(null, 97.5)) if null.size else np.nan
        region_rows.append(
            {
                "locus_class": label,
                "n_loci": len(loci),
                "observed_overlap_fraction": observed,
                "perm_mean": float(null.mean()) if null.size else np.nan,
                "perm_ci_low": ci_low,
                "perm_ci_high": ci_high,
                "perm_p_two_sided": p_value,
            }
        )
    region_df = pd.DataFrame(region_rows)
    region_df.to_csv(args.region_overlap_tsv, sep="\t", index=False)

    locus_clr: dict[str, pd.Series] = {}
    for label in ["introner_fixed_g1", "introner_polymorphic_g1"]:
        loci = introner_sets[label].copy()
        loci["nearest_clr"] = assign_nearest_clr(loci, clr_df)
        locus_clr[label] = loci["nearest_clr"]
    introns = introns.copy()
    introns["nearest_clr"] = assign_nearest_clr(introns, clr_df)
    locus_clr["non_introner_intron"] = introns["nearest_clr"]

    perlocus_rows = []
    for label, series in [
        ("introner_fixed_g1", locus_clr["introner_fixed_g1"]),
        ("introner_polymorphic_g1", locus_clr["introner_polymorphic_g1"]),
        ("non_introner_intron", locus_clr["non_introner_intron"]),
        ("genome_grid", clr_df["clr"]),
    ]:
        vals = pd.to_numeric(series, errors="coerce").dropna()
        perlocus_rows.append(
            {
                "locus_class": label,
                "n": int(len(vals)),
                "mean_clr": float(vals.mean()) if len(vals) else np.nan,
                "median_clr": float(vals.median()) if len(vals) else np.nan,
            }
        )
    perlocus_df = pd.DataFrame(perlocus_rows)

    comparisons = []
    grid_vals = clr_df["clr"].dropna()
    pairs = [
        ("introner_fixed_g1_vs_non_introner_intron", locus_clr["introner_fixed_g1"], locus_clr["non_introner_intron"]),
        ("introner_polymorphic_g1_vs_non_introner_intron", locus_clr["introner_polymorphic_g1"], locus_clr["non_introner_intron"]),
        ("introner_fixed_g1_vs_introner_polymorphic_g1", locus_clr["introner_fixed_g1"], locus_clr["introner_polymorphic_g1"]),
        ("introner_fixed_g1_vs_genome_grid", locus_clr["introner_fixed_g1"], grid_vals),
        ("introner_polymorphic_g1_vs_genome_grid", locus_clr["introner_polymorphic_g1"], grid_vals),
        ("non_introner_intron_vs_genome_grid", locus_clr["non_introner_intron"], grid_vals),
    ]
    for name, left, right in pairs:
        left_vals = pd.to_numeric(left, errors="coerce").dropna()
        right_vals = pd.to_numeric(right, errors="coerce").dropna()
        if len(left_vals) and len(right_vals):
            stat, p = mannwhitneyu(left_vals, right_vals, alternative="two-sided")
            comparisons.append((name, float(stat), float(p)))

    comp_df = pd.DataFrame(comparisons, columns=["comparison", "mannwhitney_u", "p_value"])
    perlocus_out = pd.concat([perlocus_df, comp_df], ignore_index=True, sort=False)
    perlocus_out.to_csv(args.perlocus_tsv, sep="\t", index=False)

    def overlap_fraction(label: str) -> float:
        rows = region_df[region_df["locus_class"] == label]
        return float(rows["observed_overlap_fraction"].iloc[0]) if not rows.empty else np.nan

    summary = pd.DataFrame(
        [
            ("n_introner_fixed_g1", len(introner_sets["introner_fixed_g1"])),
            ("n_introner_polymorphic_g1", len(introner_sets["introner_polymorphic_g1"])),
            ("n_non_introner_intron_loci", len(introns)),
            ("n_sweep_regions", len(regions)),
            ("n_grid_points", len(clr_df)),
            ("fixed_minus_polymorphic_overlap_fraction", overlap_fraction("introner_fixed_g1") - overlap_fraction("introner_polymorphic_g1")),
            ("fixed_minus_intron_overlap_fraction", overlap_fraction("introner_fixed_g1") - overlap_fraction("non_introner_intron")),
            ("polymorphic_minus_intron_overlap_fraction", overlap_fraction("introner_polymorphic_g1") - overlap_fraction("non_introner_intron")),
            ("permutations", args.permutations),
        ],
        columns=["metric", "value"],
    )
    summary.to_csv(args.summary_tsv, sep="\t", index=False)

    plot_enrichment(region_df, args.plot_enrichment_png, args.plot_enrichment_pdf)
    plot_clr_distributions(locus_clr, clr_df["clr"], args.plot_clr_png, args.plot_clr_pdf)

    print(f"Fixed G1 introners: {len(introner_sets['introner_fixed_g1'])} loci")
    print(f"Polymorphic G1 introners: {len(introner_sets['introner_polymorphic_g1'])} loci")
    print(f"Non-introner introns: {len(introns)} loci")
    print(f"Sweep regions: {len(regions)}")
    print(region_df.to_string(index=False))


if __name__ == "__main__":
    main()
