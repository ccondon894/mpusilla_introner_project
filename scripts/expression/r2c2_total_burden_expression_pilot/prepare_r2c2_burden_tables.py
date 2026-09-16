#!/usr/bin/env python3
"""Prepare collapsed R2C2 gene counts and CCMP1545-relative occupancy tables."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts" / "expression" / "negative_binomial_expression_pilot"))
from prepare_expression_pilot_data import median_ratio_size_factors, reference_gene_table  # noqa: E402

SAMPLES = {"834": "CCMP1545", "1614": "RCC1614", "1749": "RCC1749"}
GROUP1 = ("CCMP1545", "RCC1614")
CALLABLE = {1, 2}
PRIMARY_FAMILIES = {1, 2, 3}
CCMP_MT_CONTIG = "CCMP1545#0#scaffold_2"
CCMP_MT_START = 49808
CCMP_MT_END = 1730591
RCC_MT_CONTIG = "RCC1749#0#intronerless_contig_28"
RCC_MT_START = 25000
RCC_MT_END = 2148000
DEFAULT_QUANT = Path("/scratch2/russ/introner/splicing_fails/isoforms/sensitiveIsoforms")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def family_bin(value) -> str:
    if pd.isna(value):
        return "other"
    try:
        number = int(value)
    except (TypeError, ValueError):
        return "other"
    if number in PRIMARY_FAMILIES:
        return str(number)
    return "other"


def mating_type_genes(root: Path) -> set[str]:
    genes: set[str] = set()
    ccmp = reference_gene_table(root / "results" / "annotations" / "CCMP1545.gtf")
    genes.update(
        ccmp.loc[
            ccmp["contig"].eq(CCMP_MT_CONTIG)
            & ccmp["gene_start"].le(CCMP_MT_END)
            & ccmp["gene_end"].ge(CCMP_MT_START),
            "gene_id",
        ]
    )
    rcc = reference_gene_table(root / "results" / "annotations" / "RCC1749.gtf")
    genes.update(
        rcc.loc[
            rcc["contig"].eq(RCC_MT_CONTIG)
            & rcc["gene_start"].le(RCC_MT_END)
            & rcc["gene_end"].ge(RCC_MT_START),
            "gene_id",
        ]
    )
    return genes


def load_r2c2_gene_counts(root: Path, quant_dir: Path, eligible: set[str], mt_genes: set[str]):
    frames = []
    audits = []
    crosswalks = []
    for label, strain in SAMPLES.items():
        quant = pd.read_csv(quant_dir / f"09092025_{label}_Isoforms.filtered.clean.quant", sep="\t")
        quant = quant.loc[:, ~quant.columns.astype(str).str.startswith("Unnamed")]
        sqanti = pd.read_csv(
            root / "data" / f"sqanti3_output_{label}" / f"{label}_isoforms_classification.filtered.txt",
            sep="\t",
        )
        count_cols = [column for column in quant.columns if column.startswith(label + "_")]
        if len(count_cols) != 4:
            raise ValueError(f"{label} quant file does not have four replicate columns")
        if quant["Isoform"].duplicated().any() or set(quant["Isoform"]) != set(sqanti["isoform"]):
            raise ValueError(f"{label} isoform IDs do not match SQANTI")
        merged = quant.merge(
            sqanti[["isoform", "associated_gene", "chrom", "structural_category", *count_cols]],
            left_on="Isoform",
            right_on="isoform",
            validate="one_to_one",
            suffixes=("", "_sqanti"),
        )
        for column in count_cols:
            if not merged[column].eq(merged[column + "_sqanti"]).all():
                raise ValueError(f"{label} quant counts differ from SQANTI for {column}")
        merged["gene_id"] = merged["associated_gene"].astype(str).str.strip()
        merged["retained"] = merged["gene_id"].isin(eligible) & ~merged["gene_id"].isin(mt_genes)
        merged["strain"] = strain
        merged["isoform_total"] = merged[count_cols].sum(axis=1)
        crosswalks.append(
            merged[
                [
                    "strain",
                    "Isoform",
                    "gene_id",
                    "chrom",
                    "structural_category",
                    "retained",
                    "isoform_total",
                ]
            ]
        )
        kept = merged[merged["retained"]].copy()
        all_counts = kept.groupby("gene_id")[count_cols].sum()
        fsm = kept[kept["structural_category"].eq("full-splice_match")]
        fsm_counts = fsm.groupby("gene_id")[count_cols].sum().reindex(all_counts.index).fillna(0)
        major_idx = kept.groupby("gene_id")["isoform_total"].idxmax()
        major = kept.loc[major_idx].set_index("gene_id")[count_cols].reindex(all_counts.index).fillna(0)
        for gene_id in all_counts.index:
            for replicate in count_cols:
                frames.append(
                    {
                        "gene_id": gene_id,
                        "strain": strain,
                        "replicate": replicate,
                        "count_all": int(all_counts.loc[gene_id, replicate]),
                        "count_fsm": int(fsm_counts.loc[gene_id, replicate]),
                        "count_major": int(major.loc[gene_id, replicate]),
                    }
                )
        audits.append(
            {
                "strain": strain,
                "input_isoforms": len(quant),
                "retained_isoforms": len(kept),
                "retained_genes": int(all_counts.shape[0]),
                "fsm_isoforms": int((kept["structural_category"] == "full-splice_match").sum()),
                "input_counts": int(quant[count_cols].sum().sum()),
                "retained_counts": int(kept[count_cols].sum().sum()),
            }
        )
    counts = pd.DataFrame(frames)
    factors, n_norm = median_ratio_size_factors(
        counts.rename(columns={"count_all": "raw_count"})
    )
    counts["median_ratio_size_factor"] = counts["replicate"].map(factors)
    collapsed = counts.groupby(["gene_id", "strain"], as_index=False).agg(
        count_all=("count_all", "sum"),
        count_fsm=("count_fsm", "sum"),
        count_major=("count_major", "sum"),
        exposure=("median_ratio_size_factor", "sum"),
        n_replicates=("replicate", "nunique"),
    )
    collapsed["log_exposure"] = np.log(collapsed["exposure"])
    return collapsed, counts, pd.DataFrame(audits), pd.concat(crosswalks, ignore_index=True), n_norm


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
    rna = list(SAMPLES.values())
    mapping = mapping[mapping["sample"].isin(rna)].rename(columns={"sample": "strain"})
    presence = genotypes[genotypes["sample"].isin(rna)][
        ["ortholog_id", "sample", "presence", "family", "within_group_orthology_confidence"]
    ].rename(columns={"sample": "strain", "family": "sample_family"})
    pattern = genotypes.drop_duplicates("ortholog_id")[["ortholog_id", "group1_pattern"]]
    locus = mapping.merge(pattern, on="ortholog_id", how="left", validate="many_to_one")
    locus = locus.merge(presence, on=["ortholog_id", "strain"], how="left", validate="one_to_one")
    if not locus["presence_x"].eq(locus["presence_y"]).all():
        raise ValueError("Locus-map presence does not match the current genotype matrix")
    locus = locus.drop(columns=["presence_y"]).rename(columns={"presence_x": "presence"})
    locus["callable"] = locus["presence"].isin(CALLABLE)
    locus["present"] = locus["presence"].eq(1)
    locus["family_bin"] = np.where(locus["present"], locus["sample_family"].map(family_bin), pd.NA)
    locus["high_orthology"] = locus["within_group_orthology_confidence"].eq("high")
    locus["unique_common_gene"] = locus["ortholog_gene_mapping_status"].eq("one_exact_common_gene")
    locus["mating_type_gene"] = locus["common_gene_id"].isin(mt_genes) | locus["mating_type_gene"].eq(True)
    return locus


def gene_strain_occupancy(locus: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    counted = locus[
        locus["unique_common_gene"]
        & locus["common_gene_id"].notna()
        & ~locus["mating_type_gene"]
        & locus["high_orthology"]
        & locus["callable"]
    ].copy()
    counted["gene_id"] = counted["common_gene_id"]
    burden = counted.groupby(["gene_id", "strain"], as_index=False).agg(
        n_present=("present", "sum"),
        n_callable_loci=("ortholog_id", "nunique"),
    )
    family_counts = (
        counted[counted["present"]]
        .groupby(["gene_id", "strain", "family_bin"])
        .size()
        .unstack("family_bin", fill_value=0)
        .reset_index()
    )
    family_counts = family_counts.rename(
        columns={
            1: "n_present_family1",
            "1": "n_present_family1",
            2: "n_present_family2",
            "2": "n_present_family2",
            3: "n_present_family3",
            "3": "n_present_family3",
            "other": "n_present_other",
        }
    )
    for column in ["n_present_family1", "n_present_family2", "n_present_family3", "n_present_other"]:
        if column not in family_counts.columns:
            family_counts[column] = 0
    burden = burden.merge(
        family_counts[
            ["gene_id", "strain", "n_present_family1", "n_present_family2", "n_present_family3", "n_present_other"]
        ],
        on=["gene_id", "strain"],
        how="left",
    )
    for column in ["n_present_family1", "n_present_family2", "n_present_family3", "n_present_other"]:
        burden[column] = burden[column].fillna(0).astype(int)
    burden["n_present"] = burden["n_present"].astype(int)

    group1 = counted[counted["strain"].isin(GROUP1)]
    wide = group1.pivot_table(
        index=["ortholog_id", "gene_id"],
        columns="strain",
        values="presence",
        aggfunc="first",
    ).reset_index()
    for strain in GROUP1:
        if strain not in wide.columns:
            wide[strain] = np.nan
    wide["callable_pair"] = wide["CCMP1545"].isin(CALLABLE) & wide["RCC1614"].isin(CALLABLE)
    wide["n_extra"] = (wide["callable_pair"] & wide["CCMP1545"].eq(2) & wide["RCC1614"].eq(1)).astype(int)
    wide["n_missing"] = (wide["callable_pair"] & wide["CCMP1545"].eq(1) & wide["RCC1614"].eq(2)).astype(int)
    fam_map = counted.set_index(["ortholog_id", "strain"])["family_bin"].to_dict()
    wide["family_extra"] = [
        fam_map.get((row.ortholog_id, "RCC1614"), "other") if row.n_extra else pd.NA
        for row in wide.itertuples(index=False)
    ]
    wide["family_missing"] = [
        fam_map.get((row.ortholog_id, "CCMP1545"), "other") if row.n_missing else pd.NA
        for row in wide.itertuples(index=False)
    ]
    disc = wide.groupby("gene_id").agg(
        n_callable_pairs=("callable_pair", "sum"),
        n_extra=("n_extra", "sum"),
        n_missing=("n_missing", "sum"),
        n_uncallable_pair_loci=("callable_pair", lambda values: int((~values).sum())),
    )
    disc["n_discordant"] = disc["n_extra"] + disc["n_missing"]
    disc["delta_n"] = disc["n_extra"] - disc["n_missing"]
    for fam in ["1", "2", "3", "other"]:
        disc[f"n_extra_family{fam}"] = (
            wide.loc[wide["n_extra"].eq(1) & wide["family_extra"].eq(fam)].groupby("gene_id").size()
        )
        disc[f"n_missing_family{fam}"] = (
            wide.loc[wide["n_missing"].eq(1) & wide["family_missing"].eq(fam)].groupby("gene_id").size()
        )
        disc[f"delta_n_family{fam}"] = disc[f"n_extra_family{fam}"].fillna(0) - disc[f"n_missing_family{fam}"].fillna(0)
    conflicting = set(
        locus.loc[
            locus["ortholog_gene_mapping_status"].eq("conflicting_common_genes"),
            "native_gene_id",
        ].dropna()
    )
    mapped = counted.groupby("gene_id").agg(n_mapped_loci=("ortholog_id", "nunique"))
    disc = mapped.join(disc, how="outer").fillna(0)
    disc["conflicting_common_gene"] = disc.index.isin(conflicting)
    disc["exclude_from_paired"] = disc["n_uncallable_pair_loci"].gt(0) | disc["conflicting_common_gene"]
    return burden, disc.reset_index()


def build_pair_table(expression: pd.DataFrame, gene_flags: pd.DataFrame, mt_genes: set[str]) -> pd.DataFrame:
    group1 = expression[expression["strain"].isin(GROUP1)].copy()
    wide = group1.pivot(index="gene_id", columns="strain")
    pair = pd.DataFrame({"gene_id": wide.index})
    for name in ["count_all", "count_fsm", "count_major", "exposure", "n_replicates", "GC_content", "log_gene_span", "n_present"]:
        pair[f"{name}_CCMP1545"] = wide[name]["CCMP1545"].to_numpy()
        pair[f"{name}_RCC1614"] = wide[name]["RCC1614"].to_numpy()
    pair = pair.dropna(subset=["count_all_CCMP1545", "count_all_RCC1614"])
    pair = pair.merge(gene_flags, on="gene_id", how="left")
    int_cols = [
        "n_mapped_loci",
        "n_callable_pairs",
        "n_extra",
        "n_missing",
        "n_discordant",
        "n_uncallable_pair_loci",
        "delta_n",
    ]
    for fam in ["1", "2", "3", "other"]:
        int_cols.extend([f"n_extra_family{fam}", f"n_missing_family{fam}", f"delta_n_family{fam}"])
    for column in int_cols:
        if column in pair.columns:
            pair[column] = pair[column].fillna(0).astype(int)
        else:
            pair[column] = 0
    pair["conflicting_common_gene"] = pair["conflicting_common_gene"].eq(True)
    pair["exclude_from_paired"] = pair["exclude_from_paired"].eq(True)
    pair["mating_type_gene"] = pair["gene_id"].isin(mt_genes)
    pair["delta_GC"] = pair["GC_content_RCC1614"] - pair["GC_content_CCMP1545"]
    pair["log_exposure_ratio"] = np.log(pair["exposure_RCC1614"] / pair["exposure_CCMP1545"])
    pair["total_count"] = pair["count_all_CCMP1545"] + pair["count_all_RCC1614"]
    pair["delta_log_rate_all"] = np.log((pair["count_all_RCC1614"] + 0.5) / pair["exposure_RCC1614"]) - np.log(
        (pair["count_all_CCMP1545"] + 0.5) / pair["exposure_CCMP1545"]
    )
    pair["delta_log_rate_fsm"] = np.log((pair["count_fsm_RCC1614"] + 0.5) / pair["exposure_RCC1614"]) - np.log(
        (pair["count_fsm_CCMP1545"] + 0.5) / pair["exposure_CCMP1545"]
    )
    pair["delta_log_rate_major"] = np.log((pair["count_major_RCC1614"] + 0.5) / pair["exposure_RCC1614"]) - np.log(
        (pair["count_major_CCMP1545"] + 0.5) / pair["exposure_CCMP1545"]
    )
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
    parser.add_argument("--quant-dir", type=Path, default=DEFAULT_QUANT)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.project_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    mt_genes = mating_type_genes(root)
    ccmp_genes = set(reference_gene_table(root / "results" / "annotations" / "CCMP1545.gtf")["gene_id"])
    eligible = ccmp_genes - mt_genes
    collapsed, _replicates, audit, crosswalk, n_norm = load_r2c2_gene_counts(
        root, args.quant_dir, eligible, mt_genes
    )
    gc = pd.read_csv(root / "results" / "expression" / "glm_modeling" / "gc_content_by_gene_strain.csv")
    bounds = reference_gene_table(root / "results" / "annotations" / "CCMP1545.gtf")
    bounds["gene_span"] = bounds["gene_end"] - bounds["gene_start"] + 1
    collapsed = collapsed.merge(gc, on=["gene_id", "strain"], how="left")
    collapsed = collapsed.merge(bounds[["gene_id", "gene_span"]], on="gene_id", how="left")
    collapsed = collapsed.dropna(subset=["GC_content", "gene_span"]).copy()
    collapsed["log_gene_span"] = np.log(collapsed["gene_span"])

    locus = build_locus_table(root, mt_genes)
    burden, gene_flags = gene_strain_occupancy(locus)
    expression = collapsed.merge(burden, on=["gene_id", "strain"], how="left")
    for column in [
        "n_present",
        "n_callable_loci",
        "n_present_family1",
        "n_present_family2",
        "n_present_family3",
        "n_present_other",
    ]:
        expression[column] = expression[column].fillna(0).astype(int)
    pair = build_pair_table(expression, gene_flags, mt_genes)

    locus.to_csv(output_dir / "locus_status.tsv", sep="\t", index=False)
    burden.to_csv(output_dir / "gene_strain_occupancy.tsv", sep="\t", index=False)
    gene_flags.to_csv(output_dir / "gene_pair_flags.tsv", sep="\t", index=False)
    expression.to_csv(output_dir / "expression_gene_strain.tsv", sep="\t", index=False)
    pair.to_csv(output_dir / "gene_pair_table.tsv", sep="\t", index=False)
    audit.to_csv(output_dir / "mapping_audit.tsv", sep="\t", index=False)
    crosswalk.to_csv(output_dir / "isoform_gene_crosswalk.tsv", sep="\t", index=False)

    inputs = {
        "genotype_matrix": root / "results" / "genotyping" / "genotype_matrix.final.tsv",
        "locus_gene_map": root / "results" / "expression" / "isoform_analysis" / "turnover" / "locus_gene_map.tsv",
        "gc_content": root / "results" / "expression" / "glm_modeling" / "gc_content_by_gene_strain.csv",
    }
    (output_dir / "input_manifest.json").write_text(
        json.dumps({name: {"path": str(path), "sha256": sha256(path)} for name, path in inputs.items()}, indent=2)
        + "\n"
    )
    eligible_pairs = pair[pair["eligible_pair"]]
    summary = [
        "R2C2 total-burden expression pilot: data preparation",
        "====================================================",
        f"Mating-type genes excluded: {len(mt_genes)}",
        f"Median-ratio normalization genes: {n_norm}",
        f"Expression gene-strain rows: {len(expression)}",
        f"Expression genes: {expression['gene_id'].nunique()}",
        f"Eligible CCMP1545–RCC1614 pairs: {len(eligible_pairs)}",
        f"Pairs with nonzero Δn: {int((eligible_pairs['delta_n'] != 0).sum())}",
        f"Pairs with extra occupancy in RCC1614: {int(eligible_pairs['n_extra'].gt(0).sum())}",
        f"Pairs with missing occupancy in RCC1614: {int(eligible_pairs['n_missing'].gt(0).sum())}",
        "",
        "Primary expression is the sum of retained R2C2 isoforms. FSM-only and major-isoform",
        "sums are retained as sensitivities. Occupancy is current present count, not",
        "fixed/polymorphic class. Test B uses CCMP1545 as the strain reference.",
    ]
    (output_dir / "data_summary.txt").write_text("\n".join(summary) + "\n")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
