#!/usr/bin/env python3
"""Create global and per-contig Micro-C contact heatmaps from an mcool file."""

import argparse
import math
import re
from pathlib import Path

import cooler
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mcool", required=True, help="Input .mcool file")
    parser.add_argument("--sample", required=True, help="Sample name for plot titles")
    parser.add_argument("--global-resolution", type=int, default=5000)
    parser.add_argument("--detail-resolution", type=int, default=1000)
    parser.add_argument("--max-detail-chroms", type=int, default=12)
    parser.add_argument("--min-detail-size", type=int, default=100000)
    parser.add_argument("--outdir", required=True, help="Output directory")
    return parser.parse_args()


def safe_name(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def cooler_at_resolution(mcool, resolution):
    return cooler.Cooler(f"{mcool}::/resolutions/{resolution}")


def matrix_with_fallback(clr, region=None):
    selector = clr.matrix(balance=True)
    try:
        matrix = selector.fetch(region) if region else selector[:]
        if np.isfinite(matrix).any():
            return matrix, "balanced"
    except Exception:
        pass

    selector = clr.matrix(balance=False)
    matrix = selector.fetch(region) if region else selector[:]
    return matrix, "raw"


def positive_limits(matrix):
    values = np.asarray(matrix, dtype=float)
    values = values[np.isfinite(values) & (values > 0)]
    if values.size == 0:
        return 1e-6, 1.0

    vmin = float(np.nanpercentile(values, 5))
    vmax = float(np.nanpercentile(values, 99.5))
    if not math.isfinite(vmin) or vmin <= 0:
        vmin = float(values.min())
    if not math.isfinite(vmax) or vmax <= vmin:
        vmax = float(values.max())
    if vmax <= vmin:
        vmax = vmin * 10
    return vmin, vmax


def plot_heatmap(matrix, output_prefix, title, extent=None):
    vmin, vmax = positive_limits(matrix)
    fig, ax = plt.subplots(figsize=(8, 7))
    image = ax.imshow(
        matrix,
        cmap="magma",
        interpolation="nearest",
        norm=LogNorm(vmin=vmin, vmax=vmax),
        extent=extent,
    )
    ax.set_title(title)
    ax.set_xlabel("Genomic position")
    ax.set_ylabel("Genomic position")
    fig.colorbar(image, ax=ax, label="Contact frequency")
    fig.tight_layout()
    fig.savefig(f"{output_prefix}.png", dpi=220)
    fig.savefig(f"{output_prefix}.pdf")
    plt.close(fig)


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    detail_dir = outdir / "detail"
    detail_dir.mkdir(parents=True, exist_ok=True)

    global_clr = cooler_at_resolution(args.mcool, args.global_resolution)
    global_matrix, global_mode = matrix_with_fallback(global_clr)
    global_prefix = outdir / f"{args.sample}.global_{args.global_resolution}bp"
    plot_heatmap(
        global_matrix,
        global_prefix,
        f"{args.sample} genome-wide Micro-C ({args.global_resolution} bp, {global_mode})",
    )

    detail_clr = cooler_at_resolution(args.mcool, args.detail_resolution)
    chromsizes = detail_clr.chromsizes.sort_values(ascending=False)
    chromsizes = chromsizes[chromsizes >= args.min_detail_size]
    if args.max_detail_chroms > 0:
        chromsizes = chromsizes.head(args.max_detail_chroms)

    manifest = outdir / f"{args.sample}.heatmaps.tsv"
    with manifest.open("w") as handle:
        handle.write("sample\tchrom\tchrom_size\tresolution\tmode\tpng\tpdf\n")
        handle.write(
            f"{args.sample}\tgenome\t{int(detail_clr.chromsizes.sum())}\t"
            f"{args.global_resolution}\t{global_mode}\t"
            f"{global_prefix}.png\t{global_prefix}.pdf\n"
        )

        for chrom, size in chromsizes.items():
            matrix, mode = matrix_with_fallback(detail_clr, chrom)
            if matrix.size == 0 or not np.isfinite(matrix).any():
                continue
            output_prefix = detail_dir / f"{args.sample}.{safe_name(chrom)}.{args.detail_resolution}bp"
            plot_heatmap(
                matrix,
                output_prefix,
                f"{args.sample} {chrom} Micro-C ({args.detail_resolution} bp, {mode})",
                extent=[0, int(size), int(size), 0],
            )
            handle.write(
                f"{args.sample}\t{chrom}\t{int(size)}\t{args.detail_resolution}\t"
                f"{mode}\t{output_prefix}.png\t{output_prefix}.pdf\n"
            )


if __name__ == "__main__":
    main()
