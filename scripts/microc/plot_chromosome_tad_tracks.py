#!/usr/bin/env python3
"""Plot chromosome-scale Micro-C diagonal, insulation, introner, and gene tracks."""

import argparse
import math
import re
from pathlib import Path

import cooler
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


PRESENT = 1
ABSENT = 2


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mcool", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--genotype-matrix", required=True)
    parser.add_argument("--genes-gtf", required=True)
    parser.add_argument("--insulation", required=True)
    parser.add_argument("--boundaries", required=True)
    parser.add_argument("--chroms", required=True)
    parser.add_argument("--resolution", type=int, default=1000)
    parser.add_argument("--window", type=int, default=25000)
    parser.add_argument("--max-distance", type=int, default=250000)
    parser.add_argument("--introner-bin-bp", type=int, default=25000)
    parser.add_argument("--absent-max-length", type=int, default=200)
    parser.add_argument("--out-prefix", required=True)
    return parser.parse_args()


def safe_name(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def read_chrom_sizes(path):
    chrom_sizes = {}
    with open(path) as handle:
        for line in handle:
            chrom, size = line.rstrip("\n").split("\t")[:2]
            chrom_sizes[chrom] = int(size)
    return chrom_sizes


def chrom_number(chrom):
    match = re.search(r"scaffold_(\d+)$", chrom)
    if match:
        return match.group(1)
    match = re.search(r"#(\d+)#0$", chrom)
    if match:
        return match.group(1)
    return None


def resolve_chroms(requested, chrom_sizes):
    resolved = []
    for token in requested:
        token = token.strip()
        if not token:
            continue
        if token in chrom_sizes:
            resolved.append(token)
            continue
        normalized = token.lower().removeprefix("chr").removeprefix("scaffold_")
        matches = [
            chrom
            for chrom in chrom_sizes
            if chrom_number(chrom) == normalized
        ]
        if matches:
            resolved.append(max(matches, key=lambda chrom: chrom_sizes[chrom]))
    return resolved


def matrix_with_fallback(clr, chrom):
    try:
        matrix = clr.matrix(balance=True).fetch(chrom)
        if np.isfinite(matrix).any():
            return matrix, "balanced"
    except Exception:
        pass
    return clr.matrix(balance=False).fetch(chrom), "raw"


def diagonal_strip(matrix, resolution, max_distance):
    matrix = np.asarray(matrix, dtype=float)
    n_bins = matrix.shape[0]
    max_offset = min(max(1, max_distance // resolution), max(1, n_bins - 1))
    strip = np.full((max_offset, n_bins), np.nan, dtype=float)

    for offset in range(1, max_offset + 1):
        values = np.diag(matrix, k=offset)
        if values.size == 0:
            continue
        starts = np.arange(values.size)
        mids = starts + offset / 2.0
        columns = np.floor(mids).astype(int)
        valid = columns < n_bins
        strip[offset - 1, columns[valid]] = values[valid]
    return strip


def positive_norm(values):
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite) & (finite > 0)]
    if finite.size == 0:
        return LogNorm(vmin=1e-6, vmax=1.0)
    vmin = float(np.nanpercentile(finite, 5))
    vmax = float(np.nanpercentile(finite, 99.5))
    if not math.isfinite(vmin) or vmin <= 0:
        vmin = float(finite.min())
    if not math.isfinite(vmax) or vmax <= vmin:
        vmax = float(finite.max())
    if vmax <= vmin:
        vmax = vmin * 10
    return LogNorm(vmin=vmin, vmax=vmax)


def insulation_column(insulation, window):
    exact = f"log2_insulation_score_{window}"
    if exact in insulation.columns:
        return exact
    matches = [
        column for column in insulation.columns
        if column.startswith("log2_insulation_score")
    ]
    return matches[0] if matches else None


def read_boundaries(path, chrom_sizes):
    rows = []
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            if chrom not in chrom_sizes:
                continue
            rows.append(
                {
                    "chrom": chrom,
                    "start": int(fields[1]),
                    "end": int(fields[2]),
                }
            )
    return pd.DataFrame(rows)


def read_introner_loci(path, sample, chrom_sizes, absent_max_length):
    df = pd.read_csv(path, sep="\t")
    df = df.loc[df["sample"].astype(str) == sample].copy()
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["presence"] = df["presence"].astype(int)
    df["length"] = df["end"] - df["start"]
    df = df.loc[df["contig"].isin(chrom_sizes)].copy()
    df = df.loc[df["presence"].isin([PRESENT, ABSENT])].copy()
    df = df.loc[(df["presence"] == PRESENT) | (df["length"] <= absent_max_length)].copy()
    df["state"] = np.where(df["presence"] == PRESENT, "present", "absent")
    return df


def parse_gtf_attributes(attributes):
    parsed = {}
    for item in attributes.rstrip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if " " in item:
            key, value = item.split(" ", 1)
            parsed[key] = value.strip().strip('"')
        elif "=" in item:
            key, value = item.split("=", 1)
            parsed[key] = value.strip().strip('"')
    return parsed


def read_genes(path, chrom_sizes):
    gene_rows = []
    transcript_rows = []
    feature_rows_by_gene = {}

    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom = fields[0]
            if chrom not in chrom_sizes:
                continue
            feature = fields[2]
            start = max(0, int(fields[3]) - 1)
            end = min(chrom_sizes[chrom], int(fields[4]))
            if end <= start:
                continue
            attrs = parse_gtf_attributes(fields[8])
            gene_id = attrs.get("gene_id") or attrs.get("gene") or attrs.get("ID")
            if not gene_id:
                gene_id = attrs.get("transcript_id", f"{chrom}:{start}-{end}")

            row = {"gene_id": gene_id, "chrom": chrom, "start": start, "end": end}
            if feature == "gene":
                gene_rows.append(row)
            elif feature == "transcript":
                transcript_rows.append(row)
            elif feature in {"exon", "CDS"}:
                feature_rows_by_gene.setdefault(gene_id, []).append(row)

    if gene_rows:
        return pd.DataFrame(gene_rows).drop_duplicates("gene_id")
    if transcript_rows:
        return pd.DataFrame(transcript_rows).drop_duplicates("gene_id")

    rows = []
    for gene_id, records in feature_rows_by_gene.items():
        chroms = {record["chrom"] for record in records}
        if len(chroms) != 1:
            continue
        rows.append(
            {
                "gene_id": gene_id,
                "chrom": records[0]["chrom"],
                "start": min(record["start"] for record in records),
                "end": max(record["end"] for record in records),
            }
        )
    return pd.DataFrame(rows)


def binned_introner_counts(loci, chrom, chrom_size, bin_bp):
    edges = np.arange(0, chrom_size + bin_bp, bin_bp)
    if edges[-1] < chrom_size:
        edges = np.append(edges, chrom_size)
    centers = (edges[:-1] + edges[1:]) / 2
    chrom_loci = loci.loc[loci["contig"] == chrom].copy()
    chrom_loci["midpoint"] = ((chrom_loci["start"] + chrom_loci["end"]) / 2).astype(int)
    present_midpoints = chrom_loci.loc[chrom_loci["state"] == "present", "midpoint"]
    counts, _ = np.histogram(present_midpoints, bins=edges)
    return centers, counts


def binned_gene_counts(genes, chrom, chrom_size, bin_bp):
    edges = np.arange(0, chrom_size + bin_bp, bin_bp)
    if edges[-1] < chrom_size:
        edges = np.append(edges, chrom_size)
    centers = (edges[:-1] + edges[1:]) / 2
    chrom_genes = genes.loc[genes["chrom"] == chrom].copy()
    if chrom_genes.empty:
        return centers, np.zeros(len(centers), dtype=int)
    chrom_genes["midpoint"] = ((chrom_genes["start"] + chrom_genes["end"]) / 2).astype(int)
    counts, _ = np.histogram(chrom_genes["midpoint"], bins=edges)
    return centers, counts


def format_chr_label(chrom):
    number = chrom_number(chrom)
    return f"Chr{number}" if number else chrom


def plot_chromosome(
    clr,
    insulation,
    boundaries,
    loci,
    genes,
    chrom,
    chrom_size,
    sample,
    resolution,
    window,
    max_distance,
    introner_bin_bp,
    output_prefix,
):
    matrix, matrix_mode = matrix_with_fallback(clr, chrom)
    strip = diagonal_strip(matrix, resolution, max_distance)
    score_col = insulation_column(insulation, window)
    chrom_insulation = insulation.loc[insulation["chrom"] == chrom].copy()
    chrom_boundaries = boundaries.loc[boundaries["chrom"] == chrom].copy()
    introner_centers, introner_counts = binned_introner_counts(
        loci, chrom, chrom_size, introner_bin_bp
    )
    gene_centers, gene_counts = binned_gene_counts(
        genes, chrom, chrom_size, introner_bin_bp
    )

    fig = plt.figure(figsize=(11, 6.8))
    grid = fig.add_gridspec(4, 1, height_ratios=[2.25, 1.0, 0.85, 0.85], hspace=0.08)
    ax_heat = fig.add_subplot(grid[0])
    ax_ins = fig.add_subplot(grid[1], sharex=ax_heat)
    ax_introner = fig.add_subplot(grid[2], sharex=ax_heat)
    ax_gene = fig.add_subplot(grid[3], sharex=ax_heat)

    image = ax_heat.imshow(
        strip,
        origin="lower",
        aspect="auto",
        interpolation="nearest",
        cmap="magma",
        norm=positive_norm(strip),
        extent=[0, chrom_size, 0, max_distance],
    )
    ax_heat.set_ylabel("Micro-C\nseparation (bp)")
    ax_heat.set_title(
        f"{sample} {format_chr_label(chrom)}: Micro-C diagonal, insulation, introners, genes"
    )
    colorbar = fig.colorbar(image, ax=ax_heat, pad=0.012, fraction=0.035)
    colorbar.set_label(f"{matrix_mode} contacts")

    if score_col is not None and len(chrom_insulation):
        ax_ins.plot(
            chrom_insulation["start"],
            chrom_insulation[score_col],
            color="#293241",
            linewidth=0.9,
        )
    for _, boundary in chrom_boundaries.iterrows():
        midpoint = (int(boundary["start"]) + int(boundary["end"])) / 2
        ax_ins.axvline(midpoint, color="#d65f2d", linewidth=0.55, alpha=0.75)
        ax_heat.axvline(midpoint, color="white", linewidth=0.25, alpha=0.22)
    ax_ins.axhline(0, color="#888888", linewidth=0.5, linestyle=":")
    ax_ins.set_ylabel("Insulation")

    width = max(1, introner_bin_bp * 0.82)
    ax_introner.bar(
        introner_centers,
        introner_counts,
        width=width,
        color="#2d6cdf",
        alpha=0.78,
        label="Present introners",
    )
    ax_introner.set_ylabel("Introner\ncount")
    ax_introner.legend(frameon=False, loc="upper right")

    ax_gene.bar(
        gene_centers,
        gene_counts,
        width=width,
        color="#4c956c",
        alpha=0.78,
        label="Genes",
    )
    ax_gene.set_ylabel("Gene\ncount")
    ax_gene.set_xlabel("Genomic position (bp)")
    ax_gene.legend(frameon=False, loc="upper right")

    for ax in (ax_heat, ax_ins, ax_introner, ax_gene):
        ax.set_xlim(0, chrom_size)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    plt.setp(ax_heat.get_xticklabels(), visible=False)
    plt.setp(ax_ins.get_xticklabels(), visible=False)
    plt.setp(ax_introner.get_xticklabels(), visible=False)

    fig.tight_layout()
    fig.savefig(f"{output_prefix}.{safe_name(chrom)}.png", dpi=260)
    fig.savefig(f"{output_prefix}.{safe_name(chrom)}.pdf")
    plt.close(fig)


def main():
    args = parse_args()
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    chrom_sizes = read_chrom_sizes(args.chrom_sizes)
    requested = [chrom.strip() for chrom in args.chroms.split(",") if chrom.strip()]
    chroms = resolve_chroms(requested, chrom_sizes)
    if not chroms:
        raise ValueError(
            f"Could not resolve any requested chromosomes for {args.sample}: "
            + ",".join(requested)
        )

    clr = cooler.Cooler(f"{args.mcool}::/resolutions/{args.resolution}")
    insulation = pd.read_csv(args.insulation, sep="\t")
    boundaries = read_boundaries(args.boundaries, chrom_sizes)
    genes = read_genes(args.genes_gtf, chrom_sizes)
    loci = read_introner_loci(
        args.genotype_matrix,
        args.sample,
        chrom_sizes,
        args.absent_max_length,
    )

    manifest_rows = []
    for chrom in chroms:
        prefix = f"{out_prefix}.{safe_name(chrom)}"
        plot_chromosome(
            clr,
            insulation,
            boundaries,
            loci,
            genes,
            chrom,
            chrom_sizes[chrom],
            args.sample,
            args.resolution,
            args.window,
            args.max_distance,
            args.introner_bin_bp,
            out_prefix,
        )
        manifest_rows.append(
            {
                "sample": args.sample,
                "chrom": chrom,
                "chrom_label": format_chr_label(chrom),
                "chrom_size": chrom_sizes[chrom],
                "png": f"{prefix}.png",
                "pdf": f"{prefix}.pdf",
            }
        )

    pd.DataFrame(manifest_rows).to_csv(f"{out_prefix}.manifest.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
