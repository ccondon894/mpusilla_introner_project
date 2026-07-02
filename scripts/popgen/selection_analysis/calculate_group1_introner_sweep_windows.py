#!/usr/bin/env python3
"""Window-based selective sweep screen around Group 1 polymorphic introners."""

from __future__ import annotations

import argparse
import gzip
import math
import os
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
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
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--targets", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--window-sizes", required=True, help="Comma-separated bp sizes.")
    parser.add_argument("--min-callable-sites", type=int, default=20)
    parser.add_argument("--backgrounds-per-focal", type=int, default=50)
    parser.add_argument("--introner-windows", required=True)
    parser.add_argument("--background-windows", required=True)
    parser.add_argument("--pvalues", required=True)
    parser.add_argument("--plot-png", required=True)
    parser.add_argument("--plot-pdf", required=True)
    parser.add_argument("--focal-label", default="introner", help="Label for focal windows in summary plots.")
    return parser.parse_args()


def parse_csv(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def open_text(path: str | Path):
    path = str(path)
    if path.endswith((".gz", ".bgz")):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_gt(sample_field: str) -> int | None:
    gt = sample_field if len(sample_field) == 1 else sample_field.split(":", 1)[0]
    if gt in {".", "./.", ".|.", ""}:
        return None
    alt_count = 0
    for token in gt.replace("|", "/").split("/"):
        if token in {".", ""}:
            return None
        if token != "0":
            alt_count += 1
    return 1 if alt_count > 0 else 0


def stream_vcf_sites(vcf_path: str | Path, samples: list[str]) -> tuple[dict[str, dict[str, np.ndarray]], dict[str, int]]:
    """Return complete, non-multiallelic all-sites records keyed by contig."""
    by_contig: dict[str, dict[str, list]] = defaultdict(lambda: {"pos": [], "genotypes": [], "alt_count": []})
    stats = {
        "records_seen": 0,
        "records_used": 0,
        "skipped_missing": 0,
        "skipped_multiallelic": 0,
        "skipped_no_sample_header": 0,
    }
    sample_cols: list[int] | None = None

    with open_text(vcf_path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                vcf_samples = header[9:]
                missing = [sample for sample in samples if sample not in vcf_samples]
                if missing:
                    raise ValueError(f"VCF missing requested samples: {missing}")
                sample_cols = [9 + vcf_samples.index(sample) for sample in samples]
                continue
            if not line or line.startswith("#"):
                continue
            if sample_cols is None:
                stats["skipped_no_sample_header"] += 1
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(sample_cols):
                continue
            stats["records_seen"] += 1

            alt = fields[4]
            if "," in alt:
                stats["skipped_multiallelic"] += 1
                continue

            genotypes = []
            missing_gt = False
            for col_idx in sample_cols:
                gt = parse_gt(fields[col_idx])
                if gt is None:
                    missing_gt = True
                    break
                genotypes.append(gt)
            if missing_gt:
                stats["skipped_missing"] += 1
                continue

            contig = fields[0]
            pos = int(fields[1])
            alt_count = 0 if alt in {".", ""} else int(sum(genotypes))
            by_contig[contig]["pos"].append(pos)
            by_contig[contig]["genotypes"].append(genotypes)
            by_contig[contig]["alt_count"].append(alt_count)
            stats["records_used"] += 1

    arrays: dict[str, dict[str, np.ndarray]] = {}
    for contig, values in by_contig.items():
        order = np.argsort(np.array(values["pos"], dtype=np.int64))
        arrays[contig] = {
            "pos": np.array(values["pos"], dtype=np.int64)[order],
            "genotypes": np.array(values["genotypes"], dtype=np.int8)[order],
            "alt_count": np.array(values["alt_count"], dtype=np.int16)[order],
        }
    return arrays, stats


def tajimas_d(segregating_counts: np.ndarray, n_samples: int) -> tuple[float, float, float]:
    s = int(len(segregating_counts))
    if s == 0 or n_samples < 2:
        return math.nan, 0.0, 0.0
    a1 = sum(1 / i for i in range(1, n_samples))
    a2 = sum(1 / (i * i) for i in range(1, n_samples))
    b1 = (n_samples + 1) / (3 * (n_samples - 1))
    b2 = 2 * (n_samples * n_samples + n_samples + 3) / (9 * n_samples * (n_samples - 1))
    c1 = b1 - 1 / a1
    c2 = b2 - (n_samples + 2) / (a1 * n_samples) + a2 / (a1 * a1)
    e1 = c1 / a1
    e2 = c2 / (a1 * a1 + a2)
    pair_count = n_samples * (n_samples - 1) / 2
    pi_total = float(np.sum(segregating_counts * (n_samples - segregating_counts) / pair_count))
    theta_w_total = s / a1
    denom = math.sqrt(e1 * s + e2 * s * (s - 1))
    d_value = (pi_total - theta_w_total) / denom if denom > 0 else math.nan
    return d_value, pi_total, theta_w_total


def haplotype_homozygosity(genotypes: np.ndarray) -> tuple[float, float, int]:
    n_samples = genotypes.shape[1] if genotypes.ndim == 2 else 0
    if n_samples == 0:
        return math.nan, math.nan, 0
    haplotypes = [tuple(genotypes[:, sample_idx].tolist()) for sample_idx in range(n_samples)]
    freqs = sorted((count / n_samples for count in Counter(haplotypes).values()), reverse=True)
    h1 = float(sum(freq * freq for freq in freqs))
    if len(freqs) >= 2:
        h12 = float((freqs[0] + freqs[1]) ** 2 + sum(freq * freq for freq in freqs[2:]))
    else:
        h12 = h1
    return h1, h12, len(freqs)


def compute_window_stats(
    contig_data: dict[str, np.ndarray] | None,
    contig: str,
    start: int,
    end: int,
    n_samples: int,
    min_callable_sites: int,
) -> dict[str, object]:
    if contig_data is None:
        return empty_window_stats(contig, start, end, min_callable_sites)

    positions = contig_data["pos"]
    left = bisect_left(positions, start)
    right = bisect_right(positions, end)
    alt_counts = contig_data["alt_count"][left:right]
    callable_sites = int(len(alt_counts))
    segregating_mask = (alt_counts > 0) & (alt_counts < n_samples)
    segregating_counts = alt_counts[segregating_mask].astype(float)
    segregating_sites = int(len(segregating_counts))
    d_value, pi_total, theta_w_total = tajimas_d(segregating_counts, n_samples)
    pi_per_site = pi_total / callable_sites if callable_sites else math.nan
    theta_w_per_site = theta_w_total / callable_sites if callable_sites else math.nan

    seg_genotypes = contig_data["genotypes"][left:right][segregating_mask]
    h1, h12, n_haplotypes = haplotype_homozygosity(seg_genotypes)
    return {
        "contig": contig,
        "window_start": start,
        "window_end": end,
        "window_midpoint": (start + end) // 2,
        "callable_sites": callable_sites,
        "segregating_sites": segregating_sites,
        "pi_per_site": pi_per_site,
        "theta_w_per_site": theta_w_per_site,
        "tajimas_d": d_value,
        "h1": h1,
        "h12": h12,
        "n_haplotypes": n_haplotypes,
        "pass_min_callable": callable_sites >= min_callable_sites,
    }


def empty_window_stats(contig: str, start: int, end: int, min_callable_sites: int) -> dict[str, object]:
    return {
        "contig": contig,
        "window_start": start,
        "window_end": end,
        "window_midpoint": (start + end) // 2,
        "callable_sites": 0,
        "segregating_sites": 0,
        "pi_per_site": math.nan,
        "theta_w_per_site": math.nan,
        "tajimas_d": math.nan,
        "h1": math.nan,
        "h12": math.nan,
        "n_haplotypes": 0,
        "pass_min_callable": False,
    }


def centered_window(midpoint: int, size: int) -> tuple[int, int]:
    half = size // 2
    start = max(1, midpoint - half)
    end = start + size - 1
    return start, end


def overlaps_any(start: int, end: int, intervals: list[tuple[int, int]]) -> bool:
    return any(start <= interval_end and end >= interval_start for interval_start, interval_end in intervals)


def build_focal_windows(
    targets: pd.DataFrame,
    vcf_sites: dict[str, dict[str, np.ndarray]],
    window_sizes: list[int],
    n_samples: int,
    min_callable_sites: int,
) -> pd.DataFrame:
    rows = []
    for target in targets.itertuples(index=False):
        for window_size in window_sizes:
            start, end = centered_window(int(target.midpoint), window_size)
            stats = compute_window_stats(
                vcf_sites.get(str(target.contig)),
                str(target.contig),
                start,
                end,
                n_samples,
                min_callable_sites,
            )
            rows.append(
                {
                    "ortholog_id": target.ortholog_id,
                    "target_start": int(target.start),
                    "target_end": int(target.end),
                    "target_midpoint": int(target.midpoint),
                    "group1_present_count": int(target.group1_present_count),
                    "group1_absent_count": int(target.group1_absent_count),
                    "frequency_category": target.frequency_category,
                    "window_size": window_size,
                    **stats,
                }
            )
    return pd.DataFrame(rows)


def candidate_background_windows(
    vcf_sites: dict[str, dict[str, np.ndarray]],
    focal_windows: pd.DataFrame,
    window_size: int,
    n_samples: int,
    min_callable_sites: int,
) -> pd.DataFrame:
    focal_by_contig: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for row in focal_windows[focal_windows["window_size"] == window_size].itertuples(index=False):
        focal_by_contig[str(row.contig)].append((int(row.window_start), int(row.window_end)))

    rows = []
    for contig, contig_data in vcf_sites.items():
        positions = contig_data["pos"]
        if len(positions) == 0:
            continue
        contig_min = int(positions[0])
        contig_max = int(positions[-1])
        if contig_max - contig_min + 1 < window_size:
            continue
        for start in range(contig_min, contig_max - window_size + 2, window_size):
            end = start + window_size - 1
            if overlaps_any(start, end, focal_by_contig.get(contig, [])):
                continue
            stats = compute_window_stats(contig_data, contig, start, end, n_samples, min_callable_sites)
            if stats["pass_min_callable"]:
                rows.append({"window_size": window_size, **stats})
    return pd.DataFrame(rows)


def match_backgrounds(
    focal_windows: pd.DataFrame,
    background_candidates: pd.DataFrame,
    backgrounds_per_focal: int,
) -> pd.DataFrame:
    rows = []
    passed_focals = focal_windows[focal_windows["pass_min_callable"]].copy()
    if passed_focals.empty or background_candidates.empty:
        return pd.DataFrame()

    for focal in passed_focals.itertuples(index=False):
        same_size = background_candidates[background_candidates["window_size"] == focal.window_size]
        same_contig = same_size[same_size["contig"] == focal.contig]
        pool = same_contig if not same_contig.empty else same_size
        if pool.empty:
            continue
        pool = pool.assign(
            callable_delta=(pool["callable_sites"] - int(focal.callable_sites)).abs(),
            same_contig=pool["contig"].eq(focal.contig),
        ).sort_values(["callable_delta", "contig", "window_start"])
        for rank, bg in enumerate(pool.head(backgrounds_per_focal).itertuples(index=False), start=1):
            rows.append(
                {
                    "ortholog_id": focal.ortholog_id,
                    "target_contig": focal.contig,
                    "target_midpoint": focal.target_midpoint,
                    "window_size": focal.window_size,
                    "match_rank": rank,
                    "callable_delta": int(bg.callable_delta),
                    "same_contig": bool(bg.same_contig),
                    "contig": bg.contig,
                    "window_start": int(bg.window_start),
                    "window_end": int(bg.window_end),
                    "window_midpoint": int(bg.window_midpoint),
                    "callable_sites": int(bg.callable_sites),
                    "segregating_sites": int(bg.segregating_sites),
                    "pi_per_site": bg.pi_per_site,
                    "theta_w_per_site": bg.theta_w_per_site,
                    "tajimas_d": bg.tajimas_d,
                    "h1": bg.h1,
                    "h12": bg.h12,
                    "n_haplotypes": int(bg.n_haplotypes),
                    "pass_min_callable": bool(bg.pass_min_callable),
                }
            )
    return pd.DataFrame(rows)


def empirical_pvalue(values: pd.Series, observed: float, direction: str) -> float:
    values = pd.to_numeric(values, errors="coerce").dropna()
    if pd.isna(observed) or values.empty:
        return math.nan
    if direction == "low":
        extreme = int((values <= observed).sum())
    elif direction == "high":
        extreme = int((values >= observed).sum())
    else:
        raise ValueError(f"Unknown direction: {direction}")
    return (extreme + 1) / (len(values) + 1)


def compute_pvalues(focal_windows: pd.DataFrame, matched_backgrounds: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for focal in focal_windows.itertuples(index=False):
        if matched_backgrounds.empty:
            bg = pd.DataFrame()
        else:
            bg = matched_backgrounds[
                (matched_backgrounds["ortholog_id"] == focal.ortholog_id)
                & (matched_backgrounds["window_size"] == focal.window_size)
            ]
        p_pi = empirical_pvalue(bg["pi_per_site"], focal.pi_per_site, "low")
        p_theta = empirical_pvalue(bg["theta_w_per_site"], focal.theta_w_per_site, "low")
        p_tajima = empirical_pvalue(bg["tajimas_d"], focal.tajimas_d, "low")
        p_h1 = empirical_pvalue(bg["h1"], focal.h1, "high")
        p_h12 = empirical_pvalue(bg["h12"], focal.h12, "high")
        components = [p for p in [p_pi, p_tajima, p_h12] if pd.notna(p) and p > 0]
        composite = float(np.mean([-math.log10(p) for p in components])) if components else math.nan
        rows.append(
            {
                "ortholog_id": focal.ortholog_id,
                "contig": focal.contig,
                "target_midpoint": int(focal.target_midpoint),
                "window_size": int(focal.window_size),
                "callable_sites": int(focal.callable_sites),
                "segregating_sites": int(focal.segregating_sites),
                "pass_min_callable": bool(focal.pass_min_callable),
                "n_matched_backgrounds": int(len(bg)),
                "pi_per_site": focal.pi_per_site,
                "tajimas_d": focal.tajimas_d,
                "h12": focal.h12,
                "p_pi_low": p_pi,
                "p_theta_w_low": p_theta,
                "p_tajimas_d_low": p_tajima,
                "p_h1_high": p_h1,
                "p_h12_high": p_h12,
                "composite_sweep_score": composite,
            }
        )
    pvalues = pd.DataFrame(rows)
    if not pvalues.empty:
        pvalues = pvalues.sort_values(
            ["composite_sweep_score", "p_pi_low", "p_tajimas_d_low", "p_h12_high"],
            ascending=[False, True, True, True],
            na_position="last",
        )
    return pvalues


def plot_summary(
    focal: pd.DataFrame,
    background: pd.DataFrame,
    pvalues: pd.DataFrame,
    png: str,
    pdf: str,
    focal_label: str,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    metrics = [
        ("pi_per_site", "Nucleotide diversity per site"),
        ("tajimas_d", "Tajima's D"),
        ("h12", "H12 haplotype homozygosity"),
    ]
    for ax, (metric, label) in zip(axes.flat[:3], metrics):
        data = []
        labels = []
        for window_size in sorted(focal["window_size"].dropna().unique()):
            fvals = pd.to_numeric(
                focal[(focal["window_size"] == window_size) & focal["pass_min_callable"]][metric],
                errors="coerce",
            ).dropna()
            bvals = pd.to_numeric(background[background["window_size"] == window_size][metric], errors="coerce").dropna()
            if not fvals.empty:
                data.append(fvals)
                labels.append(f"{window_size // 1000}kb {focal_label}")
            if not bvals.empty:
                data.append(bvals.sample(min(len(bvals), 5000), random_state=42) if len(bvals) > 5000 else bvals)
                labels.append(f"{window_size // 1000}kb bg")
        if data:
            ax.boxplot(data, tick_labels=labels, showfliers=False)
            ax.tick_params(axis="x", labelrotation=45)
        ax.set_ylabel(label)

    ax = axes.flat[3]
    passed = pvalues[pvalues["pass_min_callable"] & pvalues["composite_sweep_score"].notna()].copy()
    if not passed.empty:
        ax.scatter(passed["target_midpoint"], passed["composite_sweep_score"], s=18, alpha=0.65)
        top = passed.head(8)
        for row in top.itertuples(index=False):
            ax.text(row.target_midpoint, row.composite_sweep_score, str(row.ortholog_id), fontsize=7)
    ax.set_xlabel("CCMP1545 coordinate")
    ax.set_ylabel("Composite sweep score")
    ax.set_title(f"Top candidate {focal_label} windows")
    fig.tight_layout()
    fig.savefig(png, dpi=250)
    fig.savefig(pdf)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    samples = parse_csv(args.samples)
    window_sizes = [int(value) for value in parse_csv(args.window_sizes)]
    for output in [args.introner_windows, args.background_windows, args.pvalues, args.plot_png, args.plot_pdf]:
        Path(output).parent.mkdir(parents=True, exist_ok=True)

    targets = pd.read_csv(args.targets, sep="\t")
    vcf_sites, vcf_stats = stream_vcf_sites(args.vcf, samples)
    print("VCF parse stats:")
    for key, value in vcf_stats.items():
        print(f"{key}\t{value}")

    focal_windows = build_focal_windows(targets, vcf_sites, window_sizes, len(samples), args.min_callable_sites)
    background_candidates = pd.concat(
        [
            candidate_background_windows(vcf_sites, focal_windows, window_size, len(samples), args.min_callable_sites)
            for window_size in window_sizes
        ],
        ignore_index=True,
    )
    matched_backgrounds = match_backgrounds(focal_windows, background_candidates, args.backgrounds_per_focal)
    pvalues = compute_pvalues(focal_windows, matched_backgrounds)

    focal_windows.to_csv(args.introner_windows, sep="\t", index=False)
    matched_backgrounds.to_csv(args.background_windows, sep="\t", index=False)
    pvalues.to_csv(args.pvalues, sep="\t", index=False)
    plot_summary(focal_windows, matched_backgrounds, pvalues, args.plot_png, args.plot_pdf, args.focal_label)

    print(f"Wrote {len(focal_windows)} focal windows")
    print(f"Wrote {len(matched_backgrounds)} matched background windows")
    print(f"Wrote {len(pvalues)} empirical p-value rows")


if __name__ == "__main__":
    main()
