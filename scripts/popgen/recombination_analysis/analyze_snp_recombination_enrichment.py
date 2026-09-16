#!/usr/bin/env python3
"""Test whether Group 1 polymorphic 4D SNPs occur in high-recombination regions.

Primary analysis
----------------
Compare weighted-mean 5 kb pyrho rates centered on complete-call polymorphic
fourfold-degenerate sites with rates centered on complete-call monomorphic 4D
sites using a two-sided Mann-Whitney U test. Results are reported both without
rate bounds and after applying the bounds used by the main recombination figure.

Spatial sensitivity analysis
----------------------------
In non-overlapping 5 kb windows, test the association between the fraction of
callable 4D sites that are polymorphic and log10 pyrho rate. A contig-wise
circular-shift permutation of the complete recombination-rate track preserves
the spatial autocorrelation of that track while breaking its alignment with SNP
diversity. This sensitivity test is preferred for inference because individual
sites within the same recombination segment are not independent.
"""

from __future__ import annotations

import argparse
import gzip
import math
import os
import sys
from dataclasses import dataclass
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
from scipy.stats import mannwhitneyu, rankdata, spearmanr


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
CANONICAL_BASES = {"A", "C", "G", "T"}


@dataclass(frozen=True)
class RateMap:
    starts: np.ndarray
    ends: np.ndarray
    rates: np.ndarray
    rate_run_ids: np.ndarray
    cumulative_area: np.ndarray
    map_start: int
    map_end: int


def parse_args() -> argparse.Namespace:
    outdir = PROJECT_ROOT / "analysis" / "snp_recombination_enrichment"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--vcf",
        type=Path,
        default=Path(
            "/scratch1/chris/introner-genotyping-pipeline/popgen_snp_data/"
            "mpusilla.snps.4d.notMT.group1.allsites.vcf.gz"
        ),
    )
    parser.add_argument(
        "--pyrho-dir",
        type=Path,
        default=Path("/scratch2/russ/mpusilla/pyrho/pyrho_output"),
    )
    parser.add_argument("--outdir", type=Path, default=outdir)
    parser.add_argument("--local-window-size", type=int, default=5000)
    parser.add_argument("--summary-window-size", type=int, default=5000)
    parser.add_argument("--min-callable-sites-per-window", type=int, default=20)
    parser.add_argument("--min-rate", type=float, default=1e-14)
    parser.add_argument("--max-rate", type=float, default=5e-11)
    parser.add_argument("--permutations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--accepted-filters",
        default=".,PASS",
        help="Comma-separated VCF FILTER values retained.",
    )
    parser.add_argument(
        "--color-guide",
        type=Path,
        default=PROJECT_ROOT / "master_figure_color_guide.tsv",
    )
    return parser.parse_args()


def open_text(path: Path):
    if str(path).endswith((".gz", ".bgz")):
        return gzip.open(path, "rt")
    return open(path)


def normalize_contig(value: str) -> str:
    for token in str(value).split("#"):
        if token.startswith("scaffold_"):
            return token
    return str(value)


def parse_gt(sample_field: str) -> int | None:
    """Parse the first FORMAT token despite malformed trailing reference fields."""
    gt = sample_field.split(":", 1)[0]
    if gt in {"", ".", "./.", ".|."}:
        return None
    tokens = gt.replace("|", "/").split("/")
    if any(token in {"", "."} for token in tokens):
        return None
    alleles = []
    for token in tokens:
        try:
            allele = int(token)
        except ValueError:
            return None
        if allele not in {0, 1}:
            return None
        alleles.append(allele)
    return 1 if any(alleles) else 0


def load_rate_maps(pyrho_dir: Path) -> dict[str, RateMap]:
    maps: dict[str, RateMap] = {}
    for path in sorted(pyrho_dir.glob("*.pyrho.out")):
        frame = pd.read_csv(path, sep="\t", header=None, names=["start", "end", "rate"])
        frame = frame.apply(pd.to_numeric, errors="coerce").dropna()
        frame = frame[frame["end"] > frame["start"]].sort_values(["start", "end"])
        if frame.empty:
            continue
        starts = frame["start"].to_numpy(dtype=np.int64)
        ends = frame["end"].to_numpy(dtype=np.int64)
        rates = frame["rate"].to_numpy(dtype=float)
        # Consecutive pyrho intervals often carry exactly the same rate. Keep
        # track of these runs so queries wholly inside one run can return that
        # rate directly, without subtracting two large cumulative integrals.
        rate_run_ids = np.cumsum(
            np.r_[True, rates[1:] != rates[:-1]], dtype=np.int64
        ) - 1
        # Long-double accumulation prevents cancellation when a 5 kb query in a
        # ~1e-25 segment is obtained by subtracting genome-scale integrals.
        areas = (ends - starts).astype(np.longdouble) * rates.astype(np.longdouble)
        cumulative = np.concatenate(
            [np.array([0.0], dtype=np.longdouble), np.cumsum(areas, dtype=np.longdouble)]
        )
        maps[path.name.replace(".pyrho.out", "")] = RateMap(
            starts=starts,
            ends=ends,
            rates=rates,
            rate_run_ids=rate_run_ids,
            cumulative_area=cumulative,
            map_start=int(starts[0]),
            map_end=int(ends[-1]),
        )
    if not maps:
        raise FileNotFoundError(f"No .pyrho.out files found in {pyrho_dir}")
    return maps


def integral_at(rate_map: RateMap, coordinates: np.ndarray) -> np.ndarray:
    """Integrate a piecewise-constant rate map up to each coordinate."""
    x = np.asarray(coordinates, dtype=np.int64)
    result = np.zeros(len(x), dtype=np.longdouble)
    above = x >= rate_map.map_end
    result[above] = rate_map.cumulative_area[-1]
    inside = (x > rate_map.map_start) & (x < rate_map.map_end)
    if inside.any():
        xi = x[inside]
        index = np.searchsorted(rate_map.starts, xi, side="right") - 1
        index = np.clip(index, 0, len(rate_map.starts) - 1)
        covered = np.clip(xi - rate_map.starts[index], 0, rate_map.ends[index] - rate_map.starts[index])
        result[inside] = rate_map.cumulative_area[index] + covered * rate_map.rates[index]
    return result


def weighted_rates(rate_map: RateMap, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    clipped_starts = np.maximum(np.asarray(starts, dtype=np.int64), rate_map.map_start)
    clipped_ends = np.minimum(np.asarray(ends, dtype=np.int64), rate_map.map_end)
    lengths = clipped_ends - clipped_starts
    rates = np.full(len(lengths), np.nan, dtype=float)
    valid = lengths > 0
    if valid.any():
        valid_indices = np.flatnonzero(valid)
        valid_starts = clipped_starts[valid]
        valid_ends = clipped_ends[valid]
        start_segment = np.searchsorted(rate_map.ends, valid_starts, side="right")
        end_segment = np.searchsorted(rate_map.starts, valid_ends, side="left") - 1
        same_rate_run = (
            rate_map.rate_run_ids[start_segment]
            == rate_map.rate_run_ids[end_segment]
        )
        if same_rate_run.any():
            rates[valid_indices[same_rate_run]] = rate_map.rates[
                start_segment[same_rate_run]
            ]
        crosses_segments = ~same_rate_run
        if crosses_segments.any():
            cross_starts = valid_starts[crosses_segments]
            cross_ends = valid_ends[crosses_segments]
            area = integral_at(rate_map, cross_ends) - integral_at(rate_map, cross_starts)
            rates[valid_indices[crosses_segments]] = np.asarray(
                area / (cross_ends - cross_starts), dtype=float
            )
    return rates


def parse_callable_sites(
    vcf_path: Path,
    accepted_filters: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    records: list[dict[str, object]] = []
    audit = {
        "records_seen": 0,
        "skipped_filter": 0,
        "skipped_noncanonical_ref": 0,
        "skipped_missing_or_nonbinary_gt": 0,
        "skipped_noncanonical_polymorphism": 0,
        "callable_monomorphic": 0,
        "callable_polymorphic_snp": 0,
    }
    sample_columns: list[int] | None = None

    with open_text(vcf_path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                samples = header[9:]
                missing = [sample for sample in GROUP1_SAMPLES if sample not in samples]
                if missing:
                    raise ValueError(f"VCF missing Group 1 samples: {missing}")
                sample_columns = [9 + samples.index(sample) for sample in GROUP1_SAMPLES]
                continue
            if not line or line.startswith("#"):
                continue
            if sample_columns is None:
                raise ValueError("VCF #CHROM header was not found")

            fields = line.rstrip("\n").split("\t")
            audit["records_seen"] += 1
            if len(fields) <= max(sample_columns):
                audit["skipped_missing_or_nonbinary_gt"] += 1
                continue
            if fields[6] not in accepted_filters:
                audit["skipped_filter"] += 1
                continue
            ref, alt = fields[3].upper(), fields[4].upper()
            if ref not in CANONICAL_BASES:
                audit["skipped_noncanonical_ref"] += 1
                continue
            genotypes = [parse_gt(fields[index]) for index in sample_columns]
            if any(value is None for value in genotypes):
                audit["skipped_missing_or_nonbinary_gt"] += 1
                continue
            alt_count = int(sum(genotypes))
            polymorphic = 0 < alt_count < len(GROUP1_SAMPLES)
            if polymorphic and alt not in CANONICAL_BASES:
                audit["skipped_noncanonical_polymorphism"] += 1
                continue
            class_label = "polymorphic_snp" if polymorphic else "monomorphic"
            audit[f"callable_{class_label}"] += 1
            records.append(
                {
                    "contig": fields[0],
                    "pyrho_contig": normalize_contig(fields[0]),
                    "position_1based": int(fields[1]),
                    "position_0based": int(fields[1]) - 1,
                    "ref": ref,
                    "alt": alt,
                    "alt_count": alt_count,
                    "site_class": class_label,
                }
            )
    return pd.DataFrame(records), pd.DataFrame(
        [{"audit_metric": key, "value": value} for key, value in audit.items()]
    )


def add_local_rates(
    sites: pd.DataFrame,
    maps: dict[str, RateMap],
    window_size: int,
) -> pd.DataFrame:
    sites = sites.copy()
    sites["local_window_start"] = np.nan
    sites["local_window_end"] = np.nan
    sites["local_recombination_rate"] = np.nan
    left = window_size // 2
    right = window_size - left
    for contig, index in sites.groupby("pyrho_contig").groups.items():
        rate_map = maps.get(contig)
        if rate_map is None:
            continue
        position = sites.loc[index, "position_0based"].to_numpy(dtype=np.int64)
        starts = np.maximum(0, position - left)
        ends = position + right
        sites.loc[index, "local_window_start"] = starts
        sites.loc[index, "local_window_end"] = ends
        sites.loc[index, "local_recombination_rate"] = weighted_rates(rate_map, starts, ends)
    return sites


def mannwhitney_rows(
    sites: pd.DataFrame,
    local_window_size: int,
    min_rate: float,
    max_rate: float,
) -> pd.DataFrame:
    rows = []
    for analysis_label, apply_bounds in [("unfiltered", False), ("figure_rate_bounds", True)]:
        subset = sites.dropna(subset=["local_recombination_rate"]).copy()
        if apply_bounds:
            subset = subset[
                subset["local_recombination_rate"].between(min_rate, max_rate)
            ]
        polymorphic = subset.loc[
            subset["site_class"] == "polymorphic_snp", "local_recombination_rate"
        ]
        monomorphic = subset.loc[
            subset["site_class"] == "monomorphic", "local_recombination_rate"
        ]
        test = mannwhitneyu(polymorphic, monomorphic, alternative="two-sided", method="asymptotic")
        poly_median = float(polymorphic.median())
        mono_median = float(monomorphic.median())
        rows.append(
            {
                "analysis": analysis_label,
                "local_window_size_bp": local_window_size,
                "min_rate": min_rate if apply_bounds else np.nan,
                "max_rate": max_rate if apply_bounds else np.nan,
                "polymorphic_n": len(polymorphic),
                "monomorphic_n": len(monomorphic),
                "polymorphic_median_rate": poly_median,
                "monomorphic_median_rate": mono_median,
                "median_fold_polymorphic_vs_monomorphic": poly_median / mono_median,
                "median_percent_elevation": 100.0 * (poly_median / mono_median - 1.0),
                "mannwhitney_u": float(test.statistic),
                "mannwhitney_p_two_sided": float(test.pvalue),
            }
        )
    return pd.DataFrame(rows)


def build_window_table(
    sites: pd.DataFrame,
    maps: dict[str, RateMap],
    window_size: int,
) -> pd.DataFrame:
    sites = sites.copy()
    sites["window_index"] = sites["position_0based"] // window_size
    site_counts = (
        sites.groupby(["pyrho_contig", "window_index", "site_class"])
        .size()
        .unstack(fill_value=0)
    )
    for column in ["monomorphic", "polymorphic_snp"]:
        if column not in site_counts:
            site_counts[column] = 0
    site_counts = site_counts.reset_index()
    site_counts["callable_sites"] = site_counts["monomorphic"] + site_counts["polymorphic_snp"]
    site_counts["polymorphic_fraction"] = (
        site_counts["polymorphic_snp"] / site_counts["callable_sites"]
    )

    rows = []
    for contig, rate_map in maps.items():
        first_index = rate_map.map_start // window_size
        last_index = (rate_map.map_end - 1) // window_size
        indices = np.arange(first_index, last_index + 1, dtype=np.int64)
        starts = indices * window_size
        ends = np.minimum(starts + window_size, rate_map.map_end)
        rates = weighted_rates(rate_map, starts, ends)
        rows.append(
            pd.DataFrame(
                {
                    "pyrho_contig": contig,
                    "window_index": indices,
                    "window_start": starts,
                    "window_end": ends,
                    "recombination_rate": rates,
                }
            )
        )
    windows = pd.concat(rows, ignore_index=True)
    windows = windows.merge(site_counts, how="left", on=["pyrho_contig", "window_index"])
    for column in ["monomorphic", "polymorphic_snp", "callable_sites"]:
        windows[column] = windows[column].fillna(0).astype(int)
    windows["polymorphic_fraction"] = windows["polymorphic_fraction"].astype(float)
    return windows


def pearson_correlation(left: np.ndarray, right: np.ndarray) -> float:
    left_centered = left - left.mean()
    right_centered = right - right.mean()
    denominator = math.sqrt(float(np.sum(left_centered**2) * np.sum(right_centered**2)))
    if denominator == 0:
        return math.nan
    return float(np.sum(left_centered * right_centered) / denominator)


def circular_permutation_test(
    windows: pd.DataFrame,
    min_callable: int,
    min_rate: float,
    max_rate: float,
    permutations: int,
    seed: int,
) -> tuple[pd.DataFrame, np.ndarray]:
    windows = windows.sort_values(["pyrho_contig", "window_index"]).reset_index(drop=True)
    valid = (
        (windows["callable_sites"] >= min_callable)
        & windows["polymorphic_fraction"].notna()
        & windows["recombination_rate"].between(min_rate, max_rate)
    )
    if valid.sum() < 3:
        raise ValueError("Too few valid windows for permutation analysis")

    diversity = windows.loc[valid, "polymorphic_fraction"].to_numpy(dtype=float)
    rates = windows.loc[valid, "recombination_rate"].to_numpy(dtype=float)
    observed_spearman = float(spearmanr(diversity, rates).statistic)

    # Global ranks are shifted with the full physical rate track. This retains
    # all rate-map windows, including bins without enough callable 4D sites.
    diversity_ranks = rankdata(diversity, method="average")
    full_rate_ranks = rankdata(windows["recombination_rate"].to_numpy(dtype=float), method="average")
    observed_rank_correlation = pearson_correlation(diversity_ranks, full_rate_ranks[valid.to_numpy()])

    rng = np.random.default_rng(seed)
    contig_indices = [
        np.asarray(index, dtype=int)
        for index in windows.groupby("pyrho_contig", sort=True).groups.values()
    ]
    valid_array = valid.to_numpy()
    null = np.empty(permutations, dtype=float)
    for permutation in range(permutations):
        shifted = full_rate_ranks.copy()
        for index in contig_indices:
            if len(index) > 1:
                offset = int(rng.integers(1, len(index)))
                shifted[index] = np.roll(full_rate_ranks[index], offset)
        null[permutation] = pearson_correlation(diversity_ranks, shifted[valid_array])

    null_mean = float(np.mean(null))
    p_greater = (1 + int(np.sum(null >= observed_rank_correlation))) / (permutations + 1)
    p_two_sided = (
        1
        + int(
            np.sum(
                np.abs(null - null_mean)
                >= abs(observed_rank_correlation - null_mean)
            )
        )
    ) / (permutations + 1)
    result = pd.DataFrame(
        [
            {
                "summary_window_size_bp": int(
                    (windows["window_end"] - windows["window_start"]).max()
                ),
                "min_callable_sites": min_callable,
                "min_rate": min_rate,
                "max_rate": max_rate,
                "n_valid_windows": int(valid.sum()),
                "n_contigs": int(windows.loc[valid, "pyrho_contig"].nunique()),
                "observed_spearman_rho": observed_spearman,
                "observed_global_rank_correlation": observed_rank_correlation,
                "permutations": permutations,
                "permutation_p_positive_association": p_greater,
                "permutation_p_two_sided": p_two_sided,
                "null_mean": null_mean,
                "null_sd": float(np.std(null, ddof=1)),
                "seed": seed,
            }
        ]
    )
    return result, null


def validate(
    sites: pd.DataFrame,
    windows: pd.DataFrame,
    maps: dict[str, RateMap],
) -> pd.DataFrame:
    mapped = sites.dropna(subset=["local_recombination_rate"])
    checks = [
        ("all_site_contigs_have_pyrho_maps", set(sites["pyrho_contig"]).issubset(maps)),
        ("both_site_classes_present", set(mapped["site_class"]) == {"monomorphic", "polymorphic_snp"}),
        ("all_polymorphic_alt_counts_between_1_and_10", bool(mapped.loc[mapped["site_class"] == "polymorphic_snp", "alt_count"].between(1, 10).all())),
        ("all_monomorphic_alt_counts_are_zero_or_eleven", bool(mapped.loc[mapped["site_class"] == "monomorphic", "alt_count"].isin([0, 11]).all())),
        ("all_local_rates_nonnegative", bool((mapped["local_recombination_rate"] >= 0).all())),
        ("window_counts_equal_site_counts", int(windows["callable_sites"].sum()) == len(sites)),
    ]
    return pd.DataFrame(checks, columns=["check", "pass"])


def plot_results(
    sites: pd.DataFrame,
    windows: pd.DataFrame,
    permutation_result: pd.DataFrame,
    null: np.ndarray,
    args: argparse.Namespace,
) -> None:
    colors = load_color_guide(args.color_guide)
    palette = {
        "Polymorphic SNP": get_color("Population 1 Polymorphic", colors),
        "Callable monomorphic": get_color("Population 1 Fixed", colors),
    }
    sns.set_style("whitegrid")
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))

    bounded = sites[
        sites["local_recombination_rate"].between(args.min_rate, args.max_rate)
    ].copy()
    bounded["display_class"] = bounded["site_class"].map(
        {"polymorphic_snp": "Polymorphic SNP", "monomorphic": "Callable monomorphic"}
    )
    plot_frames = []
    rng = np.random.default_rng(args.seed)
    for label, group in bounded.groupby("display_class"):
        if len(group) > 30000:
            group = group.iloc[rng.choice(len(group), 30000, replace=False)]
        plot_frames.append(group)
    site_plot = pd.concat(plot_frames, ignore_index=True)
    site_plot["log10_rate"] = np.log10(site_plot["local_recombination_rate"])
    sns.violinplot(
        data=site_plot,
        x="display_class",
        y="log10_rate",
        order=["Callable monomorphic", "Polymorphic SNP"],
        hue="display_class",
        hue_order=["Callable monomorphic", "Polymorphic SNP"],
        palette=palette,
        inner="quartile",
        cut=0,
        density_norm="width",
        legend=False,
        ax=axes[0],
    )
    axes[0].set_xlabel("")
    axes[0].set_ylabel("log10 weighted mean pyrho rate")
    axes[0].set_title("Callable 4D sites (plot subsample)")
    axes[0].tick_params(axis="x", rotation=15)

    valid = windows[
        (windows["callable_sites"] >= args.min_callable_sites_per_window)
        & windows["recombination_rate"].between(args.min_rate, args.max_rate)
    ].copy()
    valid["log10_rate"] = np.log10(valid["recombination_rate"])
    sns.regplot(
        data=valid,
        x="log10_rate",
        y="polymorphic_fraction",
        scatter_kws={"s": 10, "alpha": 0.35, "color": get_color("Population 1", colors)},
        line_kws={"color": "black", "linewidth": 1.5},
        ax=axes[1],
    )
    axes[1].set_xlabel("log10 pyrho rate (5 kb window)")
    axes[1].set_ylabel("Polymorphic fraction of callable 4D sites")
    axes[1].set_title(
        f"Window association: Spearman ρ={permutation_result.iloc[0]['observed_spearman_rho']:.3f}"
    )

    observed = float(permutation_result.iloc[0]["observed_global_rank_correlation"])
    axes[2].hist(null, bins=45, color=get_color("Ancestor", colors), edgecolor="white")
    axes[2].axvline(observed, color=get_color("Population 1", colors), linewidth=2)
    axes[2].set_xlabel("Circular-shift rank correlation")
    axes[2].set_ylabel("Permutations")
    axes[2].set_title(
        "Spatial null: "
        f"P={permutation_result.iloc[0]['permutation_p_positive_association']:.3g}"
    )

    fig.suptitle("Population 1 SNP polymorphism and local recombination", fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(args.outdir / "snp_recombination_enrichment.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    accepted_filters = {value.strip() for value in args.accepted_filters.split(",") if value.strip()}

    maps = load_rate_maps(args.pyrho_dir)
    sites, audit = parse_callable_sites(args.vcf, accepted_filters)
    sites = add_local_rates(sites, maps, args.local_window_size)
    tests = mannwhitney_rows(
        sites,
        args.local_window_size,
        args.min_rate,
        args.max_rate,
    )
    windows = build_window_table(sites, maps, args.summary_window_size)
    permutation_result, null = circular_permutation_test(
        windows,
        args.min_callable_sites_per_window,
        args.min_rate,
        args.max_rate,
        args.permutations,
        args.seed,
    )
    validation = validate(sites, windows, maps)
    if not validation["pass"].all():
        failures = validation.loc[~validation["pass"], "check"].tolist()
        raise ValueError(f"Validation failed: {failures}")

    sites.to_csv(
        args.outdir / "snp_recombination_site_rates.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    windows.to_csv(args.outdir / "snp_recombination_5kb_windows.tsv", sep="\t", index=False)
    tests.to_csv(args.outdir / "snp_recombination_mannwhitney.tsv", sep="\t", index=False)
    permutation_result.to_csv(
        args.outdir / "snp_recombination_circular_permutation.tsv", sep="\t", index=False
    )
    pd.DataFrame({"permuted_rank_correlation": null}).to_csv(
        args.outdir / "snp_recombination_circular_permutation_null.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    audit.to_csv(args.outdir / "vcf_filtering_audit.tsv", sep="\t", index=False)
    validation.to_csv(args.outdir / "validation_summary.tsv", sep="\t", index=False)
    plot_results(sites, windows, permutation_result, null, args)

    print("VCF filtering audit:")
    print(audit.to_string(index=False))
    print("\nMann-Whitney U results:")
    print(tests.to_string(index=False))
    print("\nSpatial permutation result:")
    print(permutation_result.to_string(index=False))
    print("\nValidation:")
    print(validation.to_string(index=False))


if __name__ == "__main__":
    main()
