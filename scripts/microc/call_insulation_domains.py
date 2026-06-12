#!/usr/bin/env python3
"""Calculate insulation scores and derive boundary/domain BED files."""

import argparse
from pathlib import Path

import cooler
import cooltools
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mcool", required=True, help="Input .mcool file")
    parser.add_argument("--sample", required=True, help="Sample name for plot titles")
    parser.add_argument("--resolution", type=int, default=1000)
    parser.add_argument("--window", type=int, default=25000)
    parser.add_argument("--ignore-diags", type=int, default=2)
    parser.add_argument("--min-boundary-strength", type=float, default=0.0)
    parser.add_argument(
        "--boundary-threshold-method",
        choices=("cooltools", "li", "otsu", "none"),
        default="li",
    )
    parser.add_argument("--max-plot-chroms", type=int, default=12)
    parser.add_argument("--min-plot-size", type=int, default=100000)
    parser.add_argument("--out-prefix", required=True)
    return parser.parse_args()


def calculate_insulation(clr, window, ignore_diags):
    try:
        insulation = cooltools.insulation(clr, [window], ignore_diags=ignore_diags)
        return insulation, "balanced"
    except Exception:
        insulation = cooltools.insulation(
            clr,
            [window],
            ignore_diags=ignore_diags,
            clr_weight_name=None,
        )
        return insulation, "raw"


def get_column(columns, prefix, window):
    exact = f"{prefix}_{window}"
    if exact in columns:
        return exact
    matches = [column for column in columns if column.startswith(prefix)]
    return matches[0] if matches else None


def threshold_otsu(values, bins=256):
    """Otsu threshold for positive 1D values."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan
    if np.all(values == values[0]):
        return float(values[0])

    counts, edges = np.histogram(values, bins=bins)
    centers = (edges[:-1] + edges[1:]) / 2
    weight1 = np.cumsum(counts)
    weight2 = np.cumsum(counts[::-1])[::-1]
    mean1 = np.cumsum(counts * centers) / np.maximum(weight1, 1)
    mean2 = (np.cumsum((counts * centers)[::-1]) / np.maximum(weight2[::-1], 1))[::-1]
    variance12 = weight1[:-1] * weight2[1:] * (mean1[:-1] - mean2[1:]) ** 2
    return float(centers[:-1][np.argmax(variance12)])


def threshold_li(values, tolerance=1e-6, max_iter=100):
    """Li minimum cross-entropy threshold for positive 1D values."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    values = values[values > 0]
    if len(values) == 0:
        return np.nan
    if np.all(values == values[0]):
        return float(values[0])

    threshold = float(values.mean())
    for _ in range(max_iter):
        lower = values[values <= threshold]
        upper = values[values > threshold]
        if len(lower) == 0 or len(upper) == 0:
            break
        mean_lower = lower.mean()
        mean_upper = upper.mean()
        if mean_lower <= 0 or mean_upper <= 0 or mean_lower == mean_upper:
            break
        next_threshold = (mean_lower - mean_upper) / (
            np.log(mean_lower) - np.log(mean_upper)
        )
        if abs(next_threshold - threshold) < tolerance:
            threshold = float(next_threshold)
            break
        threshold = float(next_threshold)
    return threshold


def boundary_thresholds(insulation, window):
    strength_col = get_column(insulation.columns, "boundary_strength", window)
    boundary_col = get_column(insulation.columns, "is_boundary", window)
    thresholds = {
        "none": 0.0,
        "cooltools": np.nan,
        "li": np.nan,
        "otsu": np.nan,
    }

    if not strength_col:
        return strength_col, boundary_col, thresholds

    values = insulation[strength_col].dropna().astype(float)
    values = values[values > 0]
    if len(values):
        thresholds["li"] = float(threshold_li(values.values))
        thresholds["otsu"] = float(threshold_otsu(values.values))

    if boundary_col:
        cooltools_values = insulation.loc[
            insulation[boundary_col].fillna(False).astype(bool), strength_col
        ].dropna()
        if len(cooltools_values):
            thresholds["cooltools"] = float(cooltools_values.min())

    return strength_col, boundary_col, thresholds


def select_boundaries(insulation, window, min_strength, method):
    strength_col, boundary_col, thresholds = boundary_thresholds(insulation, window)

    if strength_col:
        if method == "cooltools" and boundary_col:
            is_boundary = insulation[boundary_col].fillna(False).astype(bool)
        else:
            is_boundary = insulation[strength_col].notna()
        method_threshold = thresholds.get(method, 0.0)
        if np.isnan(method_threshold):
            method_threshold = 0.0
        threshold = max(float(min_strength), float(method_threshold))
        boundaries = insulation[
            is_boundary
            & insulation[strength_col].notna()
            & (insulation[strength_col] >= threshold)
        ].copy()
        score_values = boundaries[strength_col].fillna(0)
    elif boundary_col:
        boundaries = insulation[insulation[boundary_col].fillna(False)].copy()
        score_values = pd.Series(np.ones(len(boundaries)), index=boundaries.index)
    else:
        boundaries = insulation.iloc[0:0].copy()
        score_values = pd.Series(dtype=float)
    return boundaries, score_values, thresholds


def write_boundary_bed(boundaries, score_values, output):
    with Path(output).open("w") as handle:
        for rank, (_, row) in enumerate(boundaries.iterrows(), start=1):
            strength = float(score_values.loc[row.name]) if len(score_values) else 0.0
            score = max(0, min(1000, int(round(strength * 100))))
            handle.write(
                f"{row['chrom']}\t{int(row['start'])}\t{int(row['end'])}\t"
                f"{row['chrom']}_boundary_{rank}\t{score}\t.\t{strength:.6g}\n"
            )


def write_threshold_outputs(insulation, out_prefix, window, min_strength, selected_method):
    rows = []
    selected_boundaries = None

    for method in ("none", "cooltools", "li", "otsu"):
        boundaries, score_values, thresholds = select_boundaries(
            insulation, window, min_strength, method
        )
        suffix = "candidate_boundaries" if method == "none" else f"{method}_boundaries"
        write_boundary_bed(boundaries, score_values, f"{out_prefix}.{suffix}.bed")

        threshold = thresholds.get(method, 0.0)
        if np.isnan(threshold):
            threshold = 0.0
        rows.append(
            {
                "method": method,
                "threshold": max(float(min_strength), float(threshold)),
                "n_boundaries": len(boundaries),
                "selected": method == selected_method,
            }
        )
        if method == selected_method:
            selected_boundaries = boundaries
            write_boundary_bed(boundaries, score_values, f"{out_prefix}.boundaries.bed")

    pd.DataFrame(rows).to_csv(
        f"{out_prefix}.boundary_thresholds.tsv", sep="\t", index=False
    )
    if selected_boundaries is None:
        selected_boundaries, score_values, _ = select_boundaries(
            insulation, window, min_strength, "li"
        )
        write_boundary_bed(selected_boundaries, score_values, f"{out_prefix}.boundaries.bed")
    return selected_boundaries


def write_domains(boundaries, chromsizes, output, resolution):
    domains = []
    for chrom, chrom_size in chromsizes.items():
        starts = [0]
        ends = []
        chrom_boundaries = boundaries[boundaries["chrom"] == chrom].sort_values("start")
        for _, boundary in chrom_boundaries.iterrows():
            midpoint = int((int(boundary["start"]) + int(boundary["end"])) / 2)
            if midpoint <= starts[-1] or midpoint >= chrom_size:
                continue
            ends.append(midpoint)
            starts.append(midpoint)
        ends.append(int(chrom_size))

        for start, end in zip(starts, ends):
            if end - start < resolution:
                continue
            domains.append((chrom, int(start), int(end)))

    with Path(output).open("w") as handle:
        for idx, (chrom, start, end) in enumerate(domains, start=1):
            handle.write(f"{chrom}\t{start}\t{end}\t{chrom}_domain_{idx}\t0\t.\n")


def plot_insulation(
    insulation,
    boundaries,
    clr,
    out_prefix,
    sample,
    window,
    max_chroms,
    min_size,
    mode,
    threshold_method,
):
    score_col = get_column(insulation.columns, "log2_insulation_score", window)
    if score_col is None:
        return

    chromsizes = clr.chromsizes.sort_values(ascending=False)
    chromsizes = chromsizes[chromsizes >= min_size]
    if max_chroms > 0:
        chromsizes = chromsizes.head(max_chroms)
    if chromsizes.empty:
        return

    fig, axes = plt.subplots(len(chromsizes), 1, figsize=(9, max(2.2, 1.8 * len(chromsizes))), sharex=False)
    axes = np.atleast_1d(axes)
    for ax, (chrom, size) in zip(axes, chromsizes.items()):
        data = insulation[insulation["chrom"] == chrom]
        ax.plot(data["start"], data[score_col], color="#293241", linewidth=0.9)
        chrom_boundaries = boundaries[boundaries["chrom"] == chrom]
        for _, boundary in chrom_boundaries.iterrows():
            ax.axvline(
                (int(boundary["start"]) + int(boundary["end"])) / 2,
                color="#d65f2d",
                linewidth=0.45,
                alpha=0.5,
            )
        ax.axhline(0, color="#888888", linewidth=0.5, linestyle=":")
        ax.set_xlim(0, int(size))
        ax.set_ylabel(chrom, rotation=0, ha="right", va="center")
    axes[-1].set_xlabel("Genomic position (bp)")
    fig.suptitle(
        f"{sample} insulation score ({window} bp window, {mode}, {threshold_method} boundaries)",
        y=0.995,
    )
    fig.tight_layout()
    fig.savefig(f"{out_prefix}.insulation.png", dpi=220)
    fig.savefig(f"{out_prefix}.insulation.pdf")
    plt.close(fig)


def main():
    args = parse_args()
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    clr = cooler.Cooler(f"{args.mcool}::/resolutions/{args.resolution}")
    insulation, mode = calculate_insulation(clr, args.window, args.ignore_diags)
    insulation.insert(0, "sample", args.sample)
    insulation.insert(1, "mode", mode)

    insulation_tsv = f"{out_prefix}.insulation.tsv"
    domains_bed = f"{out_prefix}.domains.bed"

    insulation.to_csv(insulation_tsv, sep="\t", index=False)
    boundaries = write_threshold_outputs(
        insulation,
        out_prefix,
        args.window,
        args.min_boundary_strength,
        args.boundary_threshold_method,
    )
    write_domains(boundaries, clr.chromsizes, domains_bed, args.resolution)
    plot_insulation(
        insulation,
        boundaries,
        clr,
        out_prefix,
        args.sample,
        args.window,
        args.max_plot_chroms,
        args.min_plot_size,
        mode,
        args.boundary_threshold_method,
    )


if __name__ == "__main__":
    main()
