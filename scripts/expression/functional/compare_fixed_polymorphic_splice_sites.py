#!/usr/bin/env python3
"""Compare splice-site annotations of fixed and polymorphic introners."""

from __future__ import annotations

import argparse
import math
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/scratch1/chris/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
from scipy import stats
from statsmodels.stats.contingency_tables import StratifiedTable


GROUPS = {
    "group1": (
        "CCMP1545", "RCC114", "RCC1614", "RCC1698", "RCC2482", "RCC373",
        "RCC465", "RCC629", "RCC692", "RCC693", "RCC833",
    ),
    "group2": ("RCC1749", "RCC3052"),
}
RNA_SAMPLES = ("CCMP1545", "RCC1614", "RCC1749")
SAMPLE_GROUP = {sample: group for group, samples in GROUPS.items() for sample in samples}
MT_CONTIG = "CCMP1545#0#scaffold_2"
MT_START = 49_808
MT_END = 1_730_591
COLORS = {
    "none": "#999999",
    "GT_only": "#7DBFE0",
    "GC_only": "#F0A868",
    "GT-AG": "#326b77",
    "GC-AG": "#80ae9a",
    "AT-AC": "#d55e00",
    "other": "#bc272d",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--matrix", type=Path,
        default=Path("results/genotyping/genotype_matrix.final.tsv"),
    )
    parser.add_argument(
        "--boundary-dir", type=Path,
        default=Path("results/expression/splice_junctions/introner_boundary_support"),
    )
    parser.add_argument(
        "--assembly-dir", type=Path, default=Path("results/assemblies"),
    )
    parser.add_argument(
        "--outdir", type=Path,
        default=Path("results/expression/functional/fixed_vs_polymorphic_splice_sites"),
    )
    return parser.parse_args()


def bh(p_values: pd.Series) -> pd.Series:
    values = p_values.to_numpy(float)
    order = np.argsort(values)
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty(len(values))
    adjusted[order] = np.minimum(ranked, 1.0)
    return pd.Series(adjusted, index=p_values.index)


def mt_orthologs(matrix: pd.DataFrame) -> set[str]:
    ccmp = matrix[matrix["sample"].eq("CCMP1545")]
    overlap = (
        ccmp["contig"].eq(MT_CONTIG)
        & ccmp["start"].lt(MT_END)
        & ccmp["end"].gt(MT_START)
    )
    return set(ccmp.loc[overlap, "ortholog_id"])


def locus_class(row: pd.Series, group: str) -> str:
    if row[f"{group}_pattern"] == "fixed_present":
        return "fixed"
    present = int(row[f"{group}_present_count"])
    callable_count = int(row[f"{group}_callable_count"])
    n_samples = int(row[f"{group}_n_samples"])
    if callable_count == n_samples and 0 < present < n_samples:
        return "polymorphic"
    return "excluded"


def parse_annotation(value: object) -> dict[str, object]:
    text = "" if pd.isna(value) else str(value)
    gc = re.search(r"\bGC@(\d+)", text)
    gt = re.search(r"\bGT@(\d+)", text)
    donor_match = gc or gt
    donor = "GC" if gc else ("GT" if gt else "none")
    donor_offset = int(donor_match.group(1)) if donor_match else math.nan
    acceptors = list(re.finditer(r"\bAG@(\d+)", text))
    acceptor_offset = int(acceptors[-1].group(1)) if acceptors else math.nan
    has_acceptor = bool(acceptors)
    if donor == "none":
        category = "none"
    elif has_acceptor:
        category = f"{donor}-AG"
    else:
        category = f"{donor}_only"
    return {
        "annotation_category": category,
        "donor": donor,
        "donor_offset": donor_offset,
        "has_acceptor_AG": has_acceptor,
        "acceptor_offset": acceptor_offset,
    }


def modal_value(values: pd.Series) -> object:
    counts = values.value_counts(dropna=False)
    return counts.index[0]


def build_annotation_loci(matrix: pd.DataFrame) -> pd.DataFrame:
    excluded_mt = mt_orthologs(matrix)
    rows = []
    for group, samples in GROUPS.items():
        present = matrix.loc[
            matrix["sample"].isin(samples)
            & matrix["presence"].eq(1)
            & ~matrix["ortholog_id"].isin(excluded_mt)
        ].copy()
        present["locus_class"] = present.apply(locus_class, axis=1, group=group)
        present = present[present["locus_class"].isin(["fixed", "polymorphic"])]
        parsed = pd.DataFrame(present["splice_site"].map(parse_annotation).tolist())
        present = pd.concat([present.reset_index(drop=True), parsed], axis=1)
        for ortholog_id, locus in present.groupby("ortholog_id", sort=False):
            categories = locus["annotation_category"]
            modal = str(modal_value(categories))
            donor_modal = str(modal_value(locus["donor"]))
            donor_offsets = locus["donor_offset"].dropna()
            acceptor_offsets = locus["acceptor_offset"].dropna()
            first = locus.iloc[0]
            rows.append(
                {
                    "group": group,
                    "ortholog_id": ortholog_id,
                    "locus_class": first["locus_class"],
                    "family": int(first["family"]),
                    "n_present_copies": len(locus),
                    "annotation_category": modal,
                    "donor": donor_modal,
                    "recognized_donor": donor_modal in {"GT", "GC"},
                    "GC_donor": donor_modal == "GC",
                    "acceptor_AG": modal in {"GT-AG", "GC-AG"},
                    "full_pair": modal in {"GT-AG", "GC-AG"},
                    "fully_canonical_GT_AG": modal == "GT-AG",
                    "annotation_category_count": categories.nunique(),
                    "annotation_category_concordant": categories.nunique() == 1,
                    "modal_category_fraction": float(categories.eq(modal).mean()),
                    "median_donor_offset": (
                        float(donor_offsets.median()) if len(donor_offsets) else math.nan
                    ),
                    "median_acceptor_offset": (
                        float(acceptor_offsets.median()) if len(acceptor_offsets) else math.nan
                    ),
                }
            )
    return pd.DataFrame(rows)


def fisher_and_cmh(
    data: pd.DataFrame, endpoint: str, subset_name: str
) -> dict[str, object]:
    fixed = data[data["locus_class"].eq("fixed")]
    poly = data[data["locus_class"].eq("polymorphic")]
    fixed_yes = int(fixed[endpoint].sum())
    poly_yes = int(poly[endpoint].sum())
    table = np.array(
        [[poly_yes, len(poly) - poly_yes], [fixed_yes, len(fixed) - fixed_yes]]
    )
    odds, fisher_p = stats.fisher_exact(table, alternative="two-sided")

    strata = []
    for _, family in data.groupby("family"):
        family_fixed = family[family["locus_class"].eq("fixed")]
        family_poly = family[family["locus_class"].eq("polymorphic")]
        if family_fixed.empty or family_poly.empty:
            continue
        f_yes = int(family_fixed[endpoint].sum())
        p_yes = int(family_poly[endpoint].sum())
        family_table = np.array(
            [[p_yes, len(family_poly) - p_yes],
             [f_yes, len(family_fixed) - f_yes]], dtype=float
        )
        if family_table[:, 0].sum() == 0 or family_table[:, 1].sum() == 0:
            continue
        strata.append(family_table)
    if strata:
        stratified = StratifiedTable(strata, shift_zeros=True)
        cmh_odds = float(stratified.oddsratio_pooled)
        cmh_p = float(stratified.test_null_odds().pvalue)
    else:
        cmh_odds = math.nan
        cmh_p = math.nan
    return {
        "subset": subset_name,
        "endpoint": endpoint,
        "fixed_n": len(fixed),
        "polymorphic_n": len(poly),
        "fixed_yes": fixed_yes,
        "polymorphic_yes": poly_yes,
        "fixed_fraction": fixed_yes / len(fixed),
        "polymorphic_fraction": poly_yes / len(poly),
        "fisher_odds_ratio": float(odds),
        "fisher_p": float(fisher_p),
        "family_adjusted_common_odds_ratio": cmh_odds,
        "family_adjusted_cmh_p": cmh_p,
        "family_strata": len(strata),
    }


def annotation_tests(loci: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for group, group_data in loci.groupby("group", sort=False):
        subsets = [("all_loci", group_data)]
        subsets.append(
            ("category_concordant", group_data[group_data["annotation_category_concordant"]])
        )
        for subset, data in subsets:
            for endpoint in (
                "recognized_donor", "full_pair", "fully_canonical_GT_AG"
            ):
                row = fisher_and_cmh(data, endpoint, subset)
                row["group"] = group
                rows.append(row)
            donor_data = data[data["recognized_donor"]]
            row = fisher_and_cmh(donor_data, "GC_donor", subset)
            row["group"] = group
            rows.append(row)
    results = pd.DataFrame(rows)
    for subset in results["subset"].unique():
        mask = results["subset"].eq(subset)
        results.loc[mask, "fisher_fdr_bh"] = bh(results.loc[mask, "fisher_p"])
        results.loc[mask, "cmh_fdr_bh"] = bh(
            results.loc[mask, "family_adjusted_cmh_p"]
        )
    return results


COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1].upper()


def junction_motif(
    fasta: pysam.FastaFile, contig: str, start: int, end: int, strand: str
) -> tuple[str, str]:
    if strand == "+":
        return fasta.fetch(contig, start, start + 2).upper(), fasta.fetch(
            contig, end - 2, end
        ).upper()
    if strand == "-":
        return reverse_complement(fasta.fetch(contig, end - 2, end)), reverse_complement(
            fasta.fetch(contig, start, start + 2)
        )
    raise ValueError(f"Unsupported junction strand: {strand}")


def build_rna_supported_motifs(
    matrix: pd.DataFrame, boundary_dir: Path, assembly_dir: Path
) -> pd.DataFrame:
    excluded_mt = mt_orthologs(matrix)
    rows = []
    join_columns = [
        "ortholog_id", "family",
        "group1_pattern", "group1_present_count", "group1_callable_count", "group1_n_samples",
        "group2_pattern", "group2_present_count", "group2_callable_count", "group2_n_samples",
    ]
    for sample in RNA_SAMPLES:
        boundary = pd.read_csv(
            boundary_dir / f"{sample}.per_locus.tsv", sep="\t", low_memory=False
        )
        boundary = boundary.loc[
            boundary["best_junction_strand"].isin(["+", "-"])
            & boundary["max_abs_boundary_delta"].le(2)
            & ~boundary["ortholog_id"].isin(excluded_mt)
        ].copy()
        genotype = matrix.loc[matrix["sample"].eq(sample), join_columns].drop_duplicates(
            "ortholog_id"
        )
        boundary = boundary.merge(genotype, on="ortholog_id", how="left", suffixes=("", "_gt"))
        group = SAMPLE_GROUP[sample]
        boundary["locus_class"] = boundary.apply(locus_class, axis=1, group=group)
        boundary = boundary[boundary["locus_class"].isin(["fixed", "polymorphic"])]
        fasta = pysam.FastaFile(str(assembly_dir / f"{sample}.vg_paths.fa"))
        for row in boundary.itertuples():
            donor, acceptor = junction_motif(
                fasta, row.contig, int(row.best_junction_start),
                int(row.best_junction_end), row.best_junction_strand,
            )
            motif = f"{donor}-{acceptor}"
            rows.append(
                {
                    "group": group,
                    "sample": sample,
                    "ortholog_id": row.ortholog_id,
                    "locus_class": row.locus_class,
                    "family": int(row.family),
                    "contig": row.contig,
                    "junction_start": int(row.best_junction_start),
                    "junction_end": int(row.best_junction_end),
                    "junction_strand": row.best_junction_strand,
                    "best_junction_score": float(row.best_junction_score),
                    "max_abs_boundary_delta": int(row.max_abs_boundary_delta),
                    "within_2bp_replicate_count": int(row.within_2bp_replicate_count),
                    "donor_dinucleotide": donor,
                    "acceptor_dinucleotide": acceptor,
                    "observed_motif": motif,
                    "canonical_GT_AG": motif == "GT-AG",
                    "canonical_or_GC_AG": motif in {"GT-AG", "GC-AG"},
                    "GC_AG": motif == "GC-AG",
                    "noncanonical": motif not in {"GT-AG", "GC-AG"},
                }
            )
    occurrences = pd.DataFrame(rows)
    # One independent ortholog per population: retain the strongest supported
    # junction if both Group 1 RNA-seq strains carry the same ortholog.
    return (
        occurrences.sort_values(
            ["group", "ortholog_id", "best_junction_score"],
            ascending=[True, True, False],
        )
        .drop_duplicates(["group", "ortholog_id"], keep="first")
        .reset_index(drop=True)
    )


def rna_tests(motifs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for group, group_data in motifs.groupby("group", sort=False):
        for delta in (0, 1, 2):
            data = group_data[group_data["max_abs_boundary_delta"].le(delta)]
            for endpoint in ("canonical_or_GC_AG", "canonical_GT_AG"):
                row = fisher_and_cmh(data, endpoint, f"within_{delta}bp")
                row["group"] = group
                rows.append(row)
            canonical = data[data["canonical_or_GC_AG"]]
            row = fisher_and_cmh(canonical, "GC_AG", f"within_{delta}bp")
            row["group"] = group
            rows.append(row)
    results = pd.DataFrame(rows)
    results["fisher_fdr_bh"] = bh(results["fisher_p"])
    results["cmh_fdr_bh"] = bh(results["family_adjusted_cmh_p"])
    return results


def category_summary(
    data: pd.DataFrame, category: str, source: str
) -> pd.DataFrame:
    return (
        data.groupby(["group", "locus_class", category])
        .size()
        .rename("n_loci")
        .reset_index()
        .assign(source=source)
    )


def plot_categories(annotation: pd.DataFrame, rna: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    configurations = [
        (
            axes[0], annotation, "annotation_category",
            ["none", "GT_only", "GC_only", "GT-AG", "GC-AG"],
            "Sequence-search annotations",
        ),
        (
            axes[1], rna, "observed_motif",
            ["GT-AG", "GC-AG", "AT-AC", "other"],
            "RNA-supported junction motifs (within 2 bp)",
        ),
    ]
    positions = [("group1", "fixed"), ("group1", "polymorphic"),
                 ("group2", "fixed"), ("group2", "polymorphic")]
    labels = ["G1\nfixed", "G1\npoly", "G2\nfixed", "G2\npoly"]
    for ax, data, category, order, title in configurations:
        plotting = data.copy()
        if category == "observed_motif":
            plotting[category] = plotting[category].where(
                plotting[category].isin(order[:-1]), "other"
            )
        bottom = np.zeros(len(positions))
        for value in order:
            fractions = []
            for group, locus_class in positions:
                subset = plotting[
                    plotting["group"].eq(group)
                    & plotting["locus_class"].eq(locus_class)
                ]
                fractions.append(float(subset[category].eq(value).mean()))
            ax.bar(
                np.arange(len(positions)), fractions, bottom=bottom,
                color=COLORS[value], label=value.replace("_", " "), width=0.72,
            )
            bottom += np.nan_to_num(fractions)
        ax.set_title(title)
        ax.set_xticks(np.arange(len(positions)), labels)
        ax.set_ylim(0, 1)
        ax.set_ylabel("Fraction of unique ortholog loci")
        ax.legend(frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_results(
    path: Path,
    annotation: pd.DataFrame,
    annotation_stats: pd.DataFrame,
    rna: pd.DataFrame,
    rna_stats: pd.DataFrame,
) -> None:
    lines = [
        "# Fixed versus polymorphic introner splice sites",
        "",
        "Analyses use one row per ortholog locus within each population and exclude",
        "CCMP1545 mating-type-region orthologs. Population comparisons are separate.",
        "",
        "## Existing sequence-search annotations",
        "",
    ]
    counts = category_summary(annotation, "annotation_category", "annotation")
    for group in GROUPS:
        lines.extend([f"### {group}", "", "| Class | Loci | Recognized GT/GC donor | Full GT/GC-AG pair | GT-AG | GC-AG |", "|---|---:|---:|---:|---:|---:|"])
        for locus_class in ("fixed", "polymorphic"):
            data = annotation[
                annotation["group"].eq(group)
                & annotation["locus_class"].eq(locus_class)
            ]
            lines.append(
                f"| {locus_class} | {len(data)} | {data['recognized_donor'].mean():.3f} | "
                f"{data['full_pair'].mean():.3f} | {data['annotation_category'].eq('GT-AG').mean():.3f} | "
                f"{data['annotation_category'].eq('GC-AG').mean():.3f} |"
            )
        lines.append("")
    primary = annotation_stats[annotation_stats["subset"].eq("all_loci")]
    lines.extend([
        "Family-adjusted Cochran-Mantel-Haenszel tests:", "",
        "| Group | Endpoint | Common OR (poly/fixed) | P | BH FDR |",
        "|---|---|---:|---:|---:|",
    ])
    for row in primary.itertuples():
        lines.append(
            f"| {row.group} | {row.endpoint} | {row.family_adjusted_common_odds_ratio:.3f} | "
            f"{row.family_adjusted_cmh_p:.4g} | {row.cmh_fdr_bh:.4g} |"
        )
    lines.extend([
        "", "## RNA-supported observed junction motifs", "",
        "Only stranded junctions whose observed boundaries are within 2 bp of the",
        "introner-derived expected interval are included; duplicated Group 1",
        "orthologs are represented by their strongest-supported junction.", "",
        "| Group | Class | Loci | GT-AG | GC-AG | Other |",
        "|---|---|---:|---:|---:|---:|",
    ])
    for group in GROUPS:
        for locus_class in ("fixed", "polymorphic"):
            data = rna[rna["group"].eq(group) & rna["locus_class"].eq(locus_class)]
            lines.append(
                f"| {group} | {locus_class} | {len(data)} | "
                f"{int(data['observed_motif'].eq('GT-AG').sum())} | "
                f"{int(data['observed_motif'].eq('GC-AG').sum())} | "
                f"{int(data['noncanonical'].sum())} |"
            )
    rna_primary = rna_stats[rna_stats["subset"].eq("within_2bp")]
    lines.extend([
        "", "Family-adjusted tests for the within-2-bp RNA-supported set:", "",
        "| Group | Endpoint | Common OR (poly/fixed) | P | BH FDR |",
        "|---|---|---:|---:|---:|",
    ])
    for row in rna_primary.itertuples():
        lines.append(
            f"| {row.group} | {row.endpoint} | "
            f"{row.family_adjusted_common_odds_ratio:.3f} | "
            f"{row.family_adjusted_cmh_p:.4g} | {row.cmh_fdr_bh:.4g} |"
        )
    lines.extend([
        "", "## Interpretation limits", "",
        "- A missing matrix annotation means that the original search did not record a GT/GC donor; it is not proof of a non-canonical splice site.",
        "- RNA-supported motif extraction directly observes dinucleotides, but covers only loci with close, stranded short-read junction support.",
        "- Group 2 RNA-supported polymorphic loci are restricted to variants present in RCC1749 because RCC3052 lacks RNA-seq data.",
        "- Family-adjusted results are primary because introner-family frequencies differ between fixed and polymorphic loci.",
    ])
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    matrix = pd.read_csv(args.matrix, sep="\t", low_memory=False)
    annotation = build_annotation_loci(matrix)
    annotation_stats = annotation_tests(annotation)
    rna = build_rna_supported_motifs(matrix, args.boundary_dir, args.assembly_dir)
    rna_stats = rna_tests(rna)

    annotation.to_csv(args.outdir / "locus_annotation_categories.tsv", sep="\t", index=False)
    annotation_stats.to_csv(args.outdir / "annotation_class_tests.tsv", sep="\t", index=False)
    rna.to_csv(args.outdir / "rna_supported_junction_motifs.tsv", sep="\t", index=False)
    rna_stats.to_csv(args.outdir / "rna_supported_motif_tests.tsv", sep="\t", index=False)
    category_summary(annotation, "annotation_category", "annotation").to_csv(
        args.outdir / "annotation_category_summary.tsv", sep="\t", index=False
    )
    category_summary(rna, "observed_motif", "rna_supported").to_csv(
        args.outdir / "rna_motif_summary.tsv", sep="\t", index=False
    )
    plot_categories(annotation, rna, args.outdir / "splice_site_class_comparison.png")
    write_results(
        args.outdir / "RESULTS.md", annotation, annotation_stats, rna, rna_stats
    )
    print(f"Sequence-annotation loci: {len(annotation)}")
    print(f"RNA-supported unique loci: {len(rna)}")
    print(f"Wrote results to {args.outdir}")


if __name__ == "__main__":
    main()
