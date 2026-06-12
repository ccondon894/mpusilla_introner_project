#!/usr/bin/env python3
"""Estimate Micro-C contact probability as a function of genomic distance."""

import argparse
from pathlib import Path

import cooler
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
    parser.add_argument("--max-distance", type=int, default=2000000)
    parser.add_argument("--log-bins", type=int, default=80)
    parser.add_argument("--out-prefix", required=True)
    return parser.parse_args()


def fetch_sparse_matrix(clr, chrom, balanced):
    try:
        matrix = clr.matrix(balance=balanced, sparse=True).fetch(chrom).tocoo()
        data = matrix.data.astype(float)
        keep = np.isfinite(data) & (data > 0)
        if keep.any():
            matrix.data = data
            return matrix.row[keep], matrix.col[keep], data[keep], "balanced" if balanced else "raw"
    except Exception:
        if balanced:
            return fetch_sparse_matrix(clr, chrom, False)
        raise

    return np.array([], dtype=int), np.array([], dtype=int), np.array([], dtype=float), (
        "balanced" if balanced else "raw"
    )


def accumulate_by_diagonal(clr, max_distance):
    binsize = clr.binsize
    max_diag = max(1, max_distance // binsize)
    contact_sums = np.zeros(max_diag + 1, dtype=float)
    observed_pixels = np.zeros(max_diag + 1, dtype=np.int64)
    possible_pairs = np.zeros(max_diag + 1, dtype=np.int64)
    mode = "balanced"

    for chrom, size in clr.chromsizes.items():
        n_bins = int(np.ceil(size / binsize))
        usable_diag = min(max_diag, max(0, n_bins - 1))
        if usable_diag < 1:
            continue
        possible_pairs[: usable_diag + 1] += np.arange(n_bins, n_bins - usable_diag - 1, -1)

        rows, cols, data, this_mode = fetch_sparse_matrix(clr, chrom, balanced=True)
        if this_mode == "raw":
            mode = "raw"
        if data.size == 0:
            continue

        diag = cols - rows
        keep = (rows <= cols) & (diag > 0) & (diag <= max_diag)
        diag = diag[keep]
        data = data[keep]
        if data.size == 0:
            continue
        contact_sums += np.bincount(diag, weights=data, minlength=max_diag + 1)[: max_diag + 1]
        observed_pixels += np.bincount(diag, minlength=max_diag + 1)[: max_diag + 1]

    distances = np.arange(max_diag + 1) * binsize
    with np.errstate(divide="ignore", invalid="ignore"):
        contact_probability = contact_sums / possible_pairs

    table = pd.DataFrame(
        {
            "sample": clr.info.get("metadata", {}).get("sample", ""),
            "distance_bp": distances,
            "diagonal": np.arange(max_diag + 1),
            "possible_pairs": possible_pairs,
            "observed_pixels": observed_pixels,
            "contact_sum": contact_sums,
            "contact_probability": contact_probability,
            "mode": mode,
        }
    )
    table = table[(table["distance_bp"] > 0) & (table["possible_pairs"] > 0)]
    return table


def log_bin_decay(decay, log_bins):
    distances = decay["distance_bp"].to_numpy()
    min_dist = distances[distances > 0].min()
    max_dist = distances.max()
    edges = np.unique(np.rint(np.geomspace(min_dist, max_dist, log_bins + 1)).astype(int))
    if edges.size < 2:
        edges = np.array([min_dist, max_dist + 1])

    rows = []
    for start, end in zip(edges[:-1], edges[1:]):
        chunk = decay[(decay["distance_bp"] >= start) & (decay["distance_bp"] < end)]
        if chunk.empty:
            continue
        possible = chunk["possible_pairs"].sum()
        contact_sum = chunk["contact_sum"].sum()
        rows.append(
            {
                "distance_start_bp": int(start),
                "distance_end_bp": int(end),
                "distance_mid_bp": float(np.sqrt(start * end)),
                "possible_pairs": int(possible),
                "observed_pixels": int(chunk["observed_pixels"].sum()),
                "contact_sum": float(contact_sum),
                "contact_probability": float(contact_sum / possible) if possible else np.nan,
                "mode": chunk["mode"].iloc[0],
            }
        )
    return pd.DataFrame(rows)


def plot_decay(binned, out_prefix, sample, resolution):
    fig, ax = plt.subplots(figsize=(7, 5))
    plot_data = binned[np.isfinite(binned["contact_probability"]) & (binned["contact_probability"] > 0)]
    ax.plot(
        plot_data["distance_mid_bp"],
        plot_data["contact_probability"],
        marker="o",
        markersize=3,
        linewidth=1.2,
        color="#2f6f7e",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Genomic distance (bp)")
    ax.set_ylabel("Mean contact frequency")
    mode = plot_data["mode"].iloc[0] if len(plot_data) else "unknown"
    ax.set_title(f"{sample} Micro-C distance decay ({resolution} bp, {mode})")
    ax.grid(True, which="both", linewidth=0.3, alpha=0.35)
    fig.tight_layout()
    fig.savefig(f"{out_prefix}.png", dpi=220)
    fig.savefig(f"{out_prefix}.pdf")
    plt.close(fig)


def main():
    args = parse_args()
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    clr = cooler.Cooler(f"{args.mcool}::/resolutions/{args.resolution}")
    decay = accumulate_by_diagonal(clr, args.max_distance)
    decay["sample"] = args.sample
    binned = log_bin_decay(decay, args.log_bins)
    binned.insert(0, "sample", args.sample)

    decay.to_csv(f"{out_prefix}.by_diagonal.tsv", sep="\t", index=False)
    binned.to_csv(f"{out_prefix}.binned.tsv", sep="\t", index=False)
    plot_decay(binned, out_prefix, args.sample, args.resolution)


if __name__ == "__main__":
    main()
