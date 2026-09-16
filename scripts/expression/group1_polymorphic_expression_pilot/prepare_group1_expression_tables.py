#!/usr/bin/env python3
"""Build Group 1 fixed/polymorphic occupancy tables for the expression pilot.

The expression universe is featureCounts genes with exact CCMP1545 IDs after
mating-type exclusion. Introner occupancy is counted from the final genotype
matrix using group1_pattern, not CCMP1545-relative gain/loss and not RCC1749.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
SIBLING = ROOT / "scripts" / "expression" / "negative_binomial_expression_pilot"
sys.path.insert(0, str(SIBLING))

from prepare_expression_pilot_data import (  # noqa: E402
    REPLICATE_TO_STRAIN,
    median_ratio_size_factors,
    read_featurecounts,
    reference_gene_table,
)

GROUP1_RNA = ("CCMP1545", "RCC1614")
FIXED_PATTERN = "fixed_present"
POLYMORPHIC_PATTERN = "polymorphic"
CALLABLE = {1, 2}
MT_CONTIG = "CCMP1545#0#scaffold_2"
MT_START = 49808
MT_END = 1730591


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def occupancy_class(pattern: str) -> str:
    if pattern == FIXED_PATTERN:
        return "fixed"
    if pattern == POLYMORPHIC_PATTERN:
        return "polymorphic"
    return "other"


def load_expression(root: Path, mt_genes: set[str]) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    frames = []
    count_dir = root / "results" / "expression" / "counts"
    for replicate, strain in REPLICATE_TO_STRAIN.items():
        if strain not in GROUP1_RNA:
            continue
        frames.append(read_featurecounts(count_dir / f"{replicate}.counts.txt", replicate, strain))
    counts = pd.concat(frames, ignore_index=True)
    bounds = reference_gene_table(root / "results" / "annotations" / "CCMP1545.gtf")
    reference_ids = set(bounds["gene_id"]) - mt_genes
    counts = counts[counts["gene_id"].isin(reference_ids)].copy()
    size_factors, n_norm = median_ratio_size_factors(counts)
    counts["median_ratio_size_factor"] = counts["replicate"].map(size_factors)
    counts["effective_exon_kb"] = counts["effective_exon_length"] / 1000.0
    counts["exposure"] = counts["median_ratio_size_factor"] * counts["effective_exon_kb"]
    gc = pd.read_csv(root / "results" / "expression" / "glm_modeling" / "gc_content_by_gene_strain.csv")
    counts = counts.merge(gc, on=["gene_id", "strain"], how="left", validate="many_to_one")
    counts = counts.dropna(subset=["GC_content"]).copy()
    gene_span = bounds.copy()
    gene_span["gene_span"] = gene_span["gene_end"] - gene_span["gene_start"] + 1
    counts = counts.merge(gene_span[["gene_id", "gene_span"]], on="gene_id", how="left")
    collapsed = (
        counts.groupby(["gene_id", "strain"], as_index=False)
        .agg(
            raw_count=("raw_count", "sum"),
            exposure=("exposure", "sum"),
            n_replicates=("replicate", "nunique"),
            effective_exon_length=("effective_exon_length", "first"),
            GC_content=("GC_content", "first"),
            gene_span=("gene_span", "first"),
        )
    )
    collapsed["log_exposure"] = np.log(collapsed["exposure"])
    collapsed["log_gene_span"] = np.log(collapsed["gene_span"])
    collapsed["rate"] = collapsed["raw_count"] / collapsed["exposure"]
    collapsed["log_rate_pseudocount"] = np.log((collapsed["raw_count"] + 0.5) / collapsed["exposure"])
    return collapsed, counts, n_norm


def build_locus_table(root: Path, mt_genes: set[str]) -> pd.DataFrame:
    genotypes = pd.read_csv(
        root / "results" / "genotyping" / "genotype_matrix.final.tsv",
        sep="\t",
        low_memory=False,
    )
    mapping = pd.read_csv(
        root / "results" / "expression" / "isoform_analysis" / "turnover" / "locus_gene_map.tsv",
        sep="\t",
    )
    mapping = mapping[mapping["sample"].isin(GROUP1_RNA)].copy()
    mapping = mapping.rename(columns={"sample": "strain"})
    pattern = (
        genotypes.drop_duplicates("ortholog_id")
        [["ortholog_id", "group1_pattern", "within_group_orthology_confidence"]]
    )
    presence = genotypes[genotypes["sample"].isin(GROUP1_RNA)][
        ["ortholog_id", "sample", "presence"]
    ].rename(columns={"sample": "strain"})
    locus = mapping.merge(pattern, on="ortholog_id", how="left", validate="many_to_one")
    locus = locus.merge(presence, on=["ortholog_id", "strain"], how="left", validate="one_to_one")
    if not locus["presence_x"].eq(locus["presence_y"]).all():
        raise ValueError("Locus-map presence does not match the current genotype matrix")
    locus = locus.drop(columns=["presence_y"]).rename(columns={"presence_x": "presence"})
    locus["callable"] = locus["presence"].isin(CALLABLE)
    locus["present"] = locus["presence"].eq(1)
    locus["occupancy_class"] = locus["group1_pattern"].map(occupancy_class)
    locus["high_orthology"] = locus["within_group_orthology_confidence"].eq("high")
    locus["unique_common_gene"] = locus["ortholog_gene_mapping_status"].eq("one_exact_common_gene")
    locus["mating_type_gene"] = locus["common_gene_id"].isin(mt_genes) | locus["mating_type_gene"].eq(True)
    return locus


def aggregate_burden(locus: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    counted = locus[
        locus["unique_common_gene"]
        & locus["common_gene_id"].notna()
        & ~locus["mating_type_gene"]
        & locus["high_orthology"]
        & locus["callable"]
    ].copy()
    counted["gene_id"] = counted["common_gene_id"]

    keys = counted[["gene_id", "strain"]].drop_duplicates()
    burden = keys.copy()
    for occupancy, column in [
        ("fixed", "n_fixed_present"),
        ("polymorphic", "n_polymorphic_present"),
        ("other", "n_other_present"),
    ]:
        subset = counted[counted["occupancy_class"].eq(occupancy) & counted["present"]]
        counts = subset.groupby(["gene_id", "strain"], as_index=False)["ortholog_id"].nunique()
        counts = counts.rename(columns={"ortholog_id": column})
        burden = burden.merge(counts, on=["gene_id", "strain"], how="left")
        burden[column] = burden[column].fillna(0).astype(int)

    locus_n = counted.groupby(["gene_id", "strain"]).agg(
        n_callable_loci=("ortholog_id", "nunique"),
        n_polymorphic_loci=("occupancy_class", lambda values: int((values == "polymorphic").sum())),
        n_fixed_loci=("occupancy_class", lambda values: int((values == "fixed").sum())),
    )
    burden = burden.merge(locus_n.reset_index(), on=["gene_id", "strain"], how="left")
    for column in ["n_callable_loci", "n_polymorphic_loci", "n_fixed_loci"]:
        burden[column] = burden[column].fillna(0).astype(int)

    wide = counted.pivot_table(
        index=["ortholog_id", "gene_id", "occupancy_class", "group1_pattern"],
        columns="strain",
        values="presence",
        aggfunc="first",
    ).reset_index()
    for strain in GROUP1_RNA:
        if strain not in wide.columns:
            wide[strain] = np.nan
    wide["callable_pair"] = wide["CCMP1545"].isin(CALLABLE) & wide["RCC1614"].isin(CALLABLE)
    wide["discordant"] = wide["callable_pair"] & wide["CCMP1545"].ne(wide["RCC1614"])
    wide["polymorphic_discordant"] = wide["discordant"] & wide["occupancy_class"].eq("polymorphic")

    gene_flags = counted.groupby("gene_id").agg(
        n_mapped_loci=("ortholog_id", "nunique"),
        any_nonhigh_orthology=("high_orthology", lambda values: bool((~values).any())),
    )
    pair_quality = wide.groupby("gene_id").agg(
        n_callable_pairs=("callable_pair", "sum"),
        n_discordant_loci=("discordant", "sum"),
        n_polymorphic_discordant_loci=("polymorphic_discordant", "sum"),
        n_uncallable_pair_loci=("callable_pair", lambda values: int((~values).sum())),
    )
    gene_flags = gene_flags.join(pair_quality, how="outer").fillna(0)
    conflicting = set(
        locus.loc[
            locus["ortholog_gene_mapping_status"].eq("conflicting_common_genes"),
            "native_gene_id",
        ].dropna()
    )
    gene_flags["conflicting_common_gene"] = gene_flags.index.isin(conflicting)
    gene_flags["exclude_from_paired"] = (
        gene_flags["any_nonhigh_orthology"].astype(bool)
        | gene_flags["conflicting_common_gene"]
        | gene_flags["n_uncallable_pair_loci"].gt(0)
    )
    return burden, gene_flags.reset_index()


def build_pair_table(
    expression: pd.DataFrame,
    burden: pd.DataFrame,
    gene_flags: pd.DataFrame,
    mt_genes: set[str],
) -> pd.DataFrame:
    merged = expression.merge(burden, on=["gene_id", "strain"], how="left")
    count_columns = [
        "n_fixed_present",
        "n_polymorphic_present",
        "n_other_present",
        "n_callable_loci",
        "n_polymorphic_loci",
        "n_fixed_loci",
    ]
    for column in count_columns:
        merged[column] = merged[column].fillna(0).astype(int)

    wide = merged.pivot(index="gene_id", columns="strain")
    pair = pd.DataFrame({"gene_id": wide.index})
    pair = pair.set_index("gene_id")
    for name in [
        "raw_count",
        "exposure",
        "n_replicates",
        "GC_content",
        "effective_exon_length",
        "gene_span",
        *count_columns,
    ]:
        pair[f"{name}_CCMP1545"] = wide[name]["CCMP1545"]
        pair[f"{name}_RCC1614"] = wide[name]["RCC1614"]
    pair = pair.dropna(subset=["raw_count_CCMP1545", "raw_count_RCC1614"]).reset_index()
    pair = pair.merge(gene_flags, on="gene_id", how="left")
    for column in [
        "n_mapped_loci",
        "n_callable_pairs",
        "n_discordant_loci",
        "n_polymorphic_discordant_loci",
        "n_uncallable_pair_loci",
    ]:
        pair[column] = pair[column].fillna(0).astype(int)
    pair["any_nonhigh_orthology"] = pair["any_nonhigh_orthology"].eq(True)
    pair["conflicting_common_gene"] = pair["conflicting_common_gene"].eq(True)
    pair["exclude_from_paired"] = pair["exclude_from_paired"].eq(True)
    pair["delta_polymorphic_present"] = (
        pair["n_polymorphic_present_RCC1614"] - pair["n_polymorphic_present_CCMP1545"]
    )
    pair["delta_fixed_present"] = pair["n_fixed_present_RCC1614"] - pair["n_fixed_present_CCMP1545"]
    pair["delta_other_present"] = pair["n_other_present_RCC1614"] - pair["n_other_present_CCMP1545"]
    pair["delta_GC"] = pair["GC_content_RCC1614"] - pair["GC_content_CCMP1545"]
    pair["log_exposure_ratio"] = np.log(pair["exposure_RCC1614"] / pair["exposure_CCMP1545"])
    pair["total_count"] = pair["raw_count_CCMP1545"] + pair["raw_count_RCC1614"]
    pair["delta_log_rate"] = (
        np.log((pair["raw_count_RCC1614"] + 0.5) / pair["exposure_RCC1614"])
        - np.log((pair["raw_count_CCMP1545"] + 0.5) / pair["exposure_CCMP1545"])
    )
    pair["mating_type_gene"] = pair["gene_id"].isin(mt_genes)
    pair["eligible_pair"] = (
        ~pair["mating_type_gene"]
        & ~pair["exclude_from_paired"]
        & pair["n_replicates_CCMP1545"].eq(4)
        & pair["n_replicates_RCC1614"].eq(4)
        & pair["total_count"].gt(0)
        & np.isfinite(pair["log_exposure_ratio"])
    )
    return pair


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=ROOT)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.project_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    bounds = reference_gene_table(root / "results" / "annotations" / "CCMP1545.gtf")
    mt_genes = set(
        bounds.loc[
            bounds["contig"].eq(MT_CONTIG)
            & bounds["gene_start"].le(MT_END)
            & bounds["gene_end"].ge(MT_START),
            "gene_id",
        ]
    )

    expression, _replicates, n_norm = load_expression(root, mt_genes)
    locus = build_locus_table(root, mt_genes)
    burden, gene_flags = aggregate_burden(locus)
    pair = build_pair_table(expression, burden, gene_flags, mt_genes)

    merged_expression = expression.merge(burden, on=["gene_id", "strain"], how="left")
    for column in [
        "n_fixed_present",
        "n_polymorphic_present",
        "n_other_present",
        "n_callable_loci",
        "n_polymorphic_loci",
        "n_fixed_loci",
    ]:
        merged_expression[column] = merged_expression[column].fillna(0).astype(int)

    locus.to_csv(output_dir / "locus_status.tsv", sep="\t", index=False)
    burden.to_csv(output_dir / "gene_strain_burden.tsv", sep="\t", index=False)
    gene_flags.to_csv(output_dir / "gene_pair_flags.tsv", sep="\t", index=False)
    merged_expression.to_csv(output_dir / "expression_gene_strain.tsv", sep="\t", index=False)
    pair.to_csv(output_dir / "gene_pair_table.tsv", sep="\t", index=False)

    eligible = pair[pair["eligible_pair"]]
    inputs = {
        "genotype_matrix": root / "results" / "genotyping" / "genotype_matrix.final.tsv",
        "locus_gene_map": root / "results" / "expression" / "isoform_analysis" / "turnover" / "locus_gene_map.tsv",
        "gc_content": root / "results" / "expression" / "glm_modeling" / "gc_content_by_gene_strain.csv",
        "ccmp1545_gtf": root / "results" / "annotations" / "CCMP1545.gtf",
    }
    manifest = {name: {"path": str(path), "sha256": sha256(path)} for name, path in inputs.items()}
    (output_dir / "input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    summary = [
        "Group 1 polymorphic expression pilot: data preparation",
        "======================================================",
        f"Mating-type genes excluded: {len(mt_genes)}",
        f"Median-ratio normalization genes: {n_norm}",
        f"Expression gene-strain rows: {len(merged_expression)}",
        f"Expression genes: {merged_expression['gene_id'].nunique()}",
        f"Mapped unique-common-gene callable loci: {int(locus['unique_common_gene'].sum())}",
        f"Polymorphic mapped loci: {int((locus['occupancy_class'].eq('polymorphic') & locus['unique_common_gene']).sum() / 2)}",
        f"Paired genes with both strains: {len(pair)}",
        f"Eligible paired genes: {len(eligible)}",
        f"Eligible genes with a polymorphic occupancy difference: {int((eligible['delta_polymorphic_present'] != 0).sum())}",
        f"Eligible genes with a polymorphic discordant locus: {int(eligible['n_polymorphic_discordant_loci'].gt(0).sum())}",
        f"Eligible genes with a fixed occupancy difference: {int((eligible['delta_fixed_present'] != 0).sum())}",
        "",
        "Occupancy uses group1_pattern from the 11-sample genotype matrix.",
        "RCC1749 is excluded. SQANTI detection is not an inclusion criterion.",
    ]
    (output_dir / "data_summary.txt").write_text("\n".join(summary) + "\n")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
