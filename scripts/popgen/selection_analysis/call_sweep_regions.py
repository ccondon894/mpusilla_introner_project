#!/usr/bin/env python3
"""Merge SweepFinder2 CLR outputs, call sweep regions, and plot a Manhattan figure."""

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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clr-dir", required=True, help="Directory with per-contig SF2 CLR outputs.")
    parser.add_argument("--contigs-tsv", required=True, help="Contig table from build_sweepfinder_input.py.")
    parser.add_argument("--output-clr", required=True, help="Merged CLR TSV.")
    parser.add_argument("--output-regions", required=True, help="Called sweep regions BED.")
    parser.add_argument("--output-summary", required=True, help="Sweep-calling summary TSV.")
    parser.add_argument("--plot-png", required=True)
    parser.add_argument("--plot-pdf", required=True)
    parser.add_argument(
        "--site-label",
        default="4D SNPs",
        help="Site class shown in the Manhattan-plot title.",
    )
    parser.add_argument(
        "--clr-percentile",
        type=float,
        default=99.0,
        help="Empirical percentile threshold for high-CLR grid points (default top 1%%).",
    )
    parser.add_argument(
        "--region-merge-gap",
        type=int,
        default=2000,
        help="Merge adjacent high-CLR grid points separated by <= this many bp.",
    )
    parser.add_argument(
        "--introner-bed",
        default="",
        help="Optional introner BED for Manhattan rug plot.",
    )
    return parser.parse_args()


def read_sf2_clr(path: Path) -> pd.DataFrame:
    rows = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("location") or line.startswith("pos"):
                continue
            parts = line.split()
            if len(parts) < 2:
                parts = line.split("\t")
            if len(parts) < 2:
                continue
            try:
                location = float(parts[0])
                lr = float(parts[1])
            except ValueError:
                continue
            alpha = float(parts[2]) if len(parts) >= 3 else np.nan
            rows.append((int(round(location)), lr, alpha))
    if not rows:
        return pd.DataFrame(columns=["position", "clr", "alpha"])
    return pd.DataFrame(rows, columns=["position", "clr", "alpha"])


def merge_clr_outputs(clr_dir: Path, contigs: pd.DataFrame) -> pd.DataFrame:
    merged = []
    for row in contigs.itertuples(index=False):
        clr_path = clr_dir / f"{Path(row.freq_file).stem}.clr"
        if not clr_path.exists():
            stem = Path(row.grid_file).stem
            alt = clr_dir / f"{stem}.clr"
            clr_path = alt if alt.exists() else clr_path
        if not clr_path.exists():
            continue
        df = read_sf2_clr(clr_path)
        if df.empty:
            continue
        df = df.assign(contig=row.contig, contig_length=int(row.length))
        merged.append(df)
    if not merged:
        return pd.DataFrame(columns=["contig", "position", "clr", "alpha", "contig_length"])
    out = pd.concat(merged, ignore_index=True)
    return out.sort_values(["contig", "position"]).reset_index(drop=True)


def call_sweep_regions(
    clr_df: pd.DataFrame,
    percentile: float,
    merge_gap: int,
) -> tuple[pd.DataFrame, float]:
    if clr_df.empty:
        return pd.DataFrame(columns=["contig", "start", "end", "max_clr", "mean_clr", "n_grid_points"]), np.nan

    threshold = float(np.percentile(clr_df["clr"], percentile))
    high = clr_df[clr_df["clr"] >= threshold].copy()
    regions = []

    for contig, group in high.groupby("contig", sort=False):
        positions = group.sort_values("position")["position"].tolist()
        clrs = group.sort_values("position")["clr"].tolist()
        if not positions:
            continue

        start = positions[0]
        end = positions[0]
        seg_clrs = [clrs[0]]
        for idx in range(1, len(positions)):
            pos = positions[idx]
            if pos - end <= merge_gap:
                end = pos
                seg_clrs.append(clrs[idx])
            else:
                regions.append(
                    {
                        "contig": contig,
                        "start": start,
                        "end": end,
                        "max_clr": max(seg_clrs),
                        "mean_clr": float(np.mean(seg_clrs)),
                        "n_grid_points": len(seg_clrs),
                    }
                )
                start = pos
                end = pos
                seg_clrs = [clrs[idx]]
        regions.append(
            {
                "contig": contig,
                "start": start,
                "end": end,
                "max_clr": max(seg_clrs),
                "mean_clr": float(np.mean(seg_clrs)),
                "n_grid_points": len(seg_clrs),
            }
        )

    regions_df = pd.DataFrame(regions)
    if not regions_df.empty:
        regions_df = regions_df.sort_values(["contig", "start"]).reset_index(drop=True)
    return regions_df, threshold


def load_introner_rug(bed_path: str) -> pd.DataFrame:
    if not bed_path:
        return pd.DataFrame(columns=["contig", "start", "end"])
    df = pd.read_csv(bed_path, sep="\t", header=None, names=["contig", "start", "end", "name"])
    return df[["contig", "start", "end"]]


def plot_manhattan(
    clr_df: pd.DataFrame,
    regions_df: pd.DataFrame,
    threshold: float,
    introners: pd.DataFrame,
    png: str,
    pdf: str,
    site_label: str,
) -> None:
    fig, ax = plt.subplots(figsize=(14, 5))
    if clr_df.empty:
        ax.text(0.5, 0.5, "No CLR data", ha="center", va="center", transform=ax.transAxes)
        fig.savefig(png, dpi=250, bbox_inches="tight")
        fig.savefig(pdf, bbox_inches="tight")
        plt.close(fig)
        return

    contigs = clr_df["contig"].drop_duplicates().tolist()
    offsets = {}
    cumulative = 0
    for contig in contigs:
        offsets[contig] = cumulative
        cumulative += int(clr_df.loc[clr_df["contig"] == contig, "contig_length"].iloc[0]) + 50000

    x = clr_df.apply(lambda row: offsets[row["contig"]] + row["position"], axis=1)
    ax.scatter(x, clr_df["clr"], s=6, alpha=0.5, color="#2166ac", linewidths=0)

    if pd.notna(threshold):
        ax.axhline(threshold, color="#b2182b", linestyle="--", linewidth=1, label=f"CLR p{args_clr_percentile:.0f}")

    for row in regions_df.itertuples(index=False):
        x0 = offsets[row.contig] + row.start
        x1 = offsets[row.contig] + row.end
        ax.axvspan(x0, x1, color="#fdae61", alpha=0.25, linewidth=0)

    if not introners.empty:
        rug_y = ax.get_ylim()[0]
        for row in introners.itertuples(index=False):
            if row.contig not in offsets:
                continue
            midpoint = (int(row.start) + int(row.end)) // 2
            ax.plot(offsets[row.contig] + midpoint, rug_y, "|", color="#1b7837", markersize=4, alpha=0.35)

    ax.set_xlabel("Genome position (concatenated contigs)")
    ax.set_ylabel("SweepFinder2 CLR")
    ax.set_title(
        f"Group 1 SweepFinder2 scan: {site_label} (folded SFS, pyrho map)\n"
        "Note: within-G1 substructure (F_ST ~ 0.62) may inflate CLR peaks"
    )
    if pd.notna(threshold):
        ax.legend(loc="upper right", frameon=False)
    fig.tight_layout()
    fig.savefig(png, dpi=250, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)


# module-level for plot label (set in main)
args_clr_percentile = 99.0


def main() -> None:
    global args_clr_percentile
    args = parse_args()
    args_clr_percentile = args.clr_percentile

    for output in [args.output_clr, args.output_regions, args.output_summary, args.plot_png, args.plot_pdf]:
        Path(output).parent.mkdir(parents=True, exist_ok=True)

    contigs = pd.read_csv(args.contigs_tsv, sep="\t")
    clr_df = merge_clr_outputs(Path(args.clr_dir), contigs)
    regions_df, threshold = call_sweep_regions(clr_df, args.clr_percentile, args.region_merge_gap)

    clr_out = clr_df.rename(columns={"position": "position", "clr": "clr"})
    clr_out.to_csv(args.output_clr, sep="\t", index=False)

    if not regions_df.empty:
        bed = regions_df.copy()
        bed["bed_start"] = bed["start"] - 1
        bed[["contig", "bed_start", "end", "max_clr"]].to_csv(
            args.output_regions, sep="\t", header=False, index=False
        )
    else:
        Path(args.output_regions).write_text("", encoding="utf-8")

    summary = pd.DataFrame(
        [
            ("n_grid_points", len(clr_df)),
            ("n_contigs", clr_df["contig"].nunique() if not clr_df.empty else 0),
            ("clr_percentile", args.clr_percentile),
            ("clr_threshold", threshold),
            ("n_sweep_regions", len(regions_df)),
            ("region_merge_gap", args.region_merge_gap),
            ("max_clr", float(clr_df["clr"].max()) if not clr_df.empty else np.nan),
            ("median_clr", float(clr_df["clr"].median()) if not clr_df.empty else np.nan),
        ],
        columns=["metric", "value"],
    )
    summary.to_csv(args.output_summary, sep="\t", index=False)

    introners = load_introner_rug(args.introner_bed)
    plot_manhattan(
        clr_df,
        regions_df,
        threshold,
        introners,
        args.plot_png,
        args.plot_pdf,
        args.site_label,
    )

    print(f"Merged {len(clr_df)} CLR grid points across {clr_df['contig'].nunique() if not clr_df.empty else 0} contigs")
    print(f"CLR threshold (p{args.clr_percentile}): {threshold}")
    print(f"Called {len(regions_df)} sweep regions")


if __name__ == "__main__":
    main()
