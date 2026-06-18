#!/usr/bin/env python3
"""
Generate summary statistics for the final genotype matrix.

Reports ortholog group counts, within-clade insertion-site coherence,
per-clade presence patterns, cross-group ancestry (ancestral vs independent),
and supporting gene/family composition tables.
"""

from __future__ import annotations

import argparse
from collections import Counter

import pandas as pd

WITHIN_STATUS_ORDER = [
    "consistent",
    "singleton",
    "low_identity",
    "discordant",
    "uncertain",
]
WITHIN_DESCRIPTIONS = {
    "consistent": (
        "same insertion site within each clade with >=2 fingerprinted members"
    ),
    "singleton": "only one presence=1 member in the entire ortholog group",
    "low_identity": "within-clade body sequence identity below threshold",
    "discordant": "multiple distinct insertion sites within at least one clade",
    "uncertain": "insufficient fingerprint evidence for site comparison",
}
CROSS_STATUS_ORDER = [
    "ancestral",
    "likely_ancestral",
    "ancestral_low_identity",
    "independent",
    "likely_independent",
    "uncertain",
]
CROSS_ORIGIN_ORDER = [
    "ancestral",
    "independent",
    "ambiguous",
    "not_applicable",
]
PATTERN_ORDER = [
    "fixed_present",
    "present_with_missing",
    "polymorphic",
    "singleton_present",
    "singleton_absent",
    "absent",
    "no_call",
]
CROSS_REASON_ORDER = [
    "ancestral_exact",
    "ancestral_aa_match",
    "different_position",
    "compatible_diff_family",
    "exact_diff_family",
    "close_diff_aa_ctx",
    "uncertain_no_keys",
    "uncertain_single_clade",
]
ANCESTRAL_STATUSES = {
    "ancestral",
    "likely_ancestral",
    "ancestral_low_identity",
}
INDEPENDENT_STATUSES = {
    "independent",
    "likely_independent",
}


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate summary statistics for the genotype matrix"
    )
    parser.add_argument("--genotype-matrix", required=True)
    parser.add_argument("--group1-samples", nargs="+", required=True)
    parser.add_argument("--group2-samples", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def clean_cross_status(value) -> str:
    if pd.isna(value) or value == "":
        return "NA"
    return str(value)


def count_table(
    counter: Counter,
    total: int,
    order: list[str] | None = None,
    denominator: int | None = None,
) -> list[str]:
    lines = []
    denom = denominator if denominator is not None else total
    seen = set()
    items = order or sorted(counter.keys())
    for key in items:
        count = counter.get(key, 0)
        if count == 0 and key in {"discordant", "uncertain"}:
            continue
        seen.add(key)
        pct = 100 * count / denom if denom else 0
        lines.append(f"  {key:<28s}  {count:>7,d}  ({pct:5.1f}%)")
    for key, count in sorted(counter.items(), key=lambda kv: (-kv[1], str(kv[0]))):
        if key in seen:
            continue
        pct = 100 * count / denom if denom else 0
        lines.append(f"  {key:<28s}  {count:>7,d}  ({pct:5.1f}%)")
    lines.append(f"  {'TOTAL':<28s}  {total:>7,d}")
    return lines


def is_known_family(value) -> bool:
    if pd.isna(value):
        return False
    try:
        return int(value) > 0
    except (ValueError, TypeError):
        return False


def count_mixed_family_groups(df: pd.DataFrame) -> int:
    carriers = df[(df["presence"] == 1) & df["family"].map(is_known_family)].copy()
    if carriers.empty:
        return 0
    carriers["family"] = carriers["family"].astype(int)
    mixed = 0
    for _, group in carriers.groupby("ortholog_id"):
        if group["family"].nunique() > 1:
            mixed += 1
    return mixed


def gene_occupancy(loci, present_df, group_samples):
    group_present = present_df[present_df["sample"].isin(group_samples)]
    loci_genes = group_present[group_present["ortholog_id"].isin(loci)]
    has_gene = loci_genes.groupby("ortholog_id")["gene"].apply(
        lambda values: any(pd.notna(v) and v != "" for v in values)
    )
    return int(has_gene.sum()), len(has_gene)


def family_frequencies(loci, present_df, group_samples):
    group_present = present_df[
        (present_df["sample"].isin(group_samples)) & (present_df["family"] != -1)
    ]
    loci_fam = group_present[group_present["ortholog_id"].isin(loci)]
    locus_family = loci_fam.groupby("ortholog_id")["family"].agg(
        lambda values: values.mode().iloc[0] if len(values.mode()) > 0 else values.iloc[0]
    )
    return Counter(locus_family)


def main():
    args = parse_arguments()

    df = pd.read_csv(args.genotype_matrix, sep="\t")
    group1_set = set(args.group1_samples)
    group2_set = set(args.group2_samples)
    all_loci = set(df["ortholog_id"].unique())
    total_groups = len(all_loci)

    present = df[df["presence"] == 1]
    g1_loci = set(present[present["sample"].isin(group1_set)]["ortholog_id"])
    g2_loci = set(present[present["sample"].isin(group2_set)]["ortholog_id"])
    shared_loci = g1_loci & g2_loci
    neither_loci = all_loci - g1_loci - g2_loci

    fixed_shared = set()
    for oid in shared_loci:
        rows = df[df["ortholog_id"] == oid]
        g1_rows = rows[rows["sample"].isin(group1_set)]
        g2_rows = rows[rows["sample"].isin(group2_set)]
        if (g1_rows["presence"] == 1).all() and (g2_rows["presence"] == 1).all():
            fixed_shared.add(oid)

    group_status = df.drop_duplicates("ortholog_id").set_index("ortholog_id")
    within_status_counts = Counter(group_status["within_group_status"].fillna(""))
    cross_status_counts = Counter(
        group_status["cross_group_status"].map(clean_cross_status)
    )

    has_derived = {
        "group1_pattern",
        "group2_pattern",
        "group1_callability",
        "group2_callability",
        "within_group_orthology_confidence",
        "cross_group_origin",
        "cross_group_confidence",
        "present_in_both_groups",
    }.issubset(group_status.columns)

    g1_pattern_counts = Counter()
    g2_pattern_counts = Counter()
    cross_origin_counts = Counter()
    cross_confidence_counts = Counter()
    within_confidence_counts = Counter()
    cross_reason_counts = Counter()
    biology_counts = Counter()

    if has_derived:
        g1_pattern_counts = Counter(group_status["group1_pattern"].fillna(""))
        g2_pattern_counts = Counter(group_status["group2_pattern"].fillna(""))
        cross_origin_counts = Counter(group_status["cross_group_origin"].fillna(""))
        cross_confidence_counts = Counter(
            group_status["cross_group_confidence"].fillna("")
        )
        within_confidence_counts = Counter(
            group_status["within_group_orthology_confidence"].fillna("")
        )

    if "cross_group_reason" in group_status.columns:
        reason_mask = group_status["cross_group_reason"].fillna("").astype(str) != ""
        cross_reason_counts = Counter(group_status.loc[reason_mask, "cross_group_reason"])

    if has_derived:
        for oid, row in group_status.iterrows():
            g1_pat = row.get("group1_pattern", "")
            g2_pat = row.get("group2_pattern", "")
            origin = row.get("cross_group_origin", "not_applicable")
            if not row.get("present_in_both_groups", False):
                biology_counts["within_clade_only"] += 1
                continue
            if origin == "ancestral":
                if g1_pat == "fixed_present" and g2_pat == "fixed_present":
                    biology_counts["shared_ancestral_fixed_both_clades"] += 1
                else:
                    biology_counts["shared_ancestral_not_fixed_both_clades"] += 1
            elif origin == "independent":
                biology_counts["shared_independent_insertion"] += 1
            elif origin == "ambiguous":
                biology_counts["shared_ambiguous_cross_group"] += 1
            else:
                biology_counts["shared_other_cross_group"] += 1

    na_groups = set(
        group_status[group_status["cross_group_status"].map(clean_cross_status) == "NA"].index
    )
    na_g1_only = len({oid for oid in na_groups if oid in g1_loci and oid not in g2_loci})
    na_g2_only = len({oid for oid in na_groups if oid in g2_loci and oid not in g1_loci})
    na_both_presence = len(
        {oid for oid in na_groups if oid in g1_loci and oid in g2_loci}
    )
    na_neither = len(
        {oid for oid in na_groups if oid not in g1_loci and oid not in g2_loci}
    )

    g1_in_gene, g1_gene_total = gene_occupancy(g1_loci, present, group1_set)
    g2_in_gene, g2_gene_total = gene_occupancy(g2_loci, present, group2_set)
    g1_fam_counts = family_frequencies(g1_loci, present, group1_set)
    g2_fam_counts = family_frequencies(g2_loci, present, group2_set)
    mixed_family_groups = count_mixed_family_groups(df)

    cross_group_total = total_groups - cross_status_counts.get("NA", 0)
    ancestral_total = sum(cross_status_counts.get(s, 0) for s in ANCESTRAL_STATUSES)
    independent_total = sum(cross_status_counts.get(s, 0) for s in INDEPENDENT_STATUSES)

    with open(args.output, "w") as handle:
        handle.write("=" * 72 + "\n")
        handle.write("GENOTYPE MATRIX SUMMARY STATISTICS\n")
        handle.write("=" * 72 + "\n\n")

        handle.write("TERMINOLOGY (KEY DISTINCTIONS)\n")
        handle.write("-" * 40 + "\n")
        handle.write(
            "within_group_status = insertion-site coherence within each clade.\n"
            "  'consistent' means carriers in G1 (and/or G2) share the same insertion\n"
            "  site by fingerprint/body-sequence evidence. It is NOT fixation.\n"
            "group1_pattern / group2_pattern = per-clade presence frequency after coverage.\n"
            "  'fixed_present' means all callable samples in that clade are present.\n"
            "cross_group_status / cross_group_origin = ancestry between clades when both\n"
            "  clades have carriers: shared ancestral insertion vs independent insertion.\n\n"
        )

        handle.write("1. ORTHOLOG GROUP COUNTS\n")
        handle.write("-" * 40 + "\n")
        handle.write(f"Total ortholog groups:        {total_groups:,}\n")
        handle.write(f"Group 1 loci (>=1 present):   {len(g1_loci):,}\n")
        handle.write(f"Group 2 loci (>=1 present):   {len(g2_loci):,}\n")
        handle.write(f"Present in both clades:       {len(shared_loci):,}\n")
        handle.write(f"Present in neither clade:     {len(neither_loci):,}\n")
        handle.write(f"Fixed present in all samples: {len(fixed_shared):,}\n")
        handle.write(f"Mixed-family carrier groups:  {mixed_family_groups:,}\n\n")

        handle.write("2. WITHIN-CLADE INSERTION-SITE CLASSIFICATION\n")
        handle.write("-" * 40 + "\n")
        handle.write("within_group_status (one row per ortholog group):\n")
        for line in count_table(within_status_counts, total_groups, WITHIN_STATUS_ORDER):
            handle.write(line + "\n")
        handle.write("\nDefinitions:\n")
        for status in WITHIN_STATUS_ORDER:
            if status in within_status_counts:
                handle.write(f"  {status}: {WITHIN_DESCRIPTIONS[status]}\n")
        handle.write("\n")

        if has_derived:
            handle.write("3. PER-CLADE PRESENCE PATTERNS (POST-COVERAGE)\n")
            handle.write("-" * 40 + "\n")
            handle.write("Group 1 pattern:\n")
            for line in count_table(g1_pattern_counts, total_groups, PATTERN_ORDER):
                handle.write(line + "\n")
            handle.write("\nGroup 2 pattern:\n")
            for line in count_table(g2_pattern_counts, total_groups, PATTERN_ORDER):
                handle.write(line + "\n")
            handle.write("\n")

            handle.write("4. WITHIN-GROUP ORTHOLOGY CONFIDENCE\n")
            handle.write("-" * 40 + "\n")
            for line in count_table(
                within_confidence_counts, total_groups, ["high", "low_identity", "discordant", "uncertain", "unknown"]
            ):
                handle.write(line + "\n")
            handle.write("\n")

        section = 5 if has_derived else 3
        handle.write(f"{section}. CROSS-GROUP ANCESTRY (G1 vs G2)\n")
        handle.write("-" * 40 + "\n")
        handle.write(
            f"Groups with both clade carriers (n={cross_group_total:,}):\n"
        )
        for line in count_table(
            cross_status_counts,
            total_groups,
            CROSS_STATUS_ORDER,
            denominator=cross_group_total if cross_group_total else total_groups,
        ):
            if line.strip().startswith("TOTAL"):
                continue
            if line.split()[0] == "NA":
                continue
            handle.write(line + "\n")
        handle.write(
            f"\nCollapsed ancestry (groups with both clade carriers):\n"
        )
        handle.write(
            f"  {'ancestral (all tiers)':<28s}  {ancestral_total:>7,d}  "
            f"({100 * ancestral_total / cross_group_total if cross_group_total else 0:5.1f}%)\n"
        )
        handle.write(
            f"  {'independent (all tiers)':<28s}  {independent_total:>7,d}  "
            f"({100 * independent_total / cross_group_total if cross_group_total else 0:5.1f}%)\n"
        )
        handle.write(
            f"  {'uncertain':<28s}  {cross_status_counts.get('uncertain', 0):>7,d}  "
            f"({100 * cross_status_counts.get('uncertain', 0) / cross_group_total if cross_group_total else 0:5.1f}%)\n"
        )
        handle.write(
            f"\nWithin-clade-only groups (cross_group_status = NA): "
            f"{cross_status_counts.get('NA', 0):,}\n"
        )
        handle.write(f"  Present only in G1: {na_g1_only:,}\n")
        handle.write(f"  Present only in G2: {na_g2_only:,}\n")
        if na_both_presence:
            handle.write(f"  Present in both*: {na_both_presence:,}\n")
            handle.write(
                "    * retained NA because cross-group comparison lacked evidence\n"
            )
        if na_neither:
            handle.write(f"  Present in neither: {na_neither:,}\n")
        handle.write("\n")

        if has_derived:
            section += 1
            handle.write(f"{section}. CROSS-GROUP ORIGIN AND CONFIDENCE\n")
            handle.write("-" * 40 + "\n")
            handle.write("cross_group_origin:\n")
            for line in count_table(cross_origin_counts, total_groups, CROSS_ORIGIN_ORDER):
                handle.write(line + "\n")
            handle.write("\ncross_group_confidence:\n")
            for line in count_table(
                cross_confidence_counts,
                total_groups,
                ["high", "moderate", "low_identity", "low", "not_applicable"],
            ):
                handle.write(line + "\n")
            handle.write("\n")

        if cross_reason_counts:
            section += 1
            handle.write(f"{section}. CROSS-GROUP REASON (FINE-GRAINED)\n")
            handle.write("-" * 40 + "\n")
            for line in count_table(
                cross_reason_counts,
                sum(cross_reason_counts.values()),
                CROSS_REASON_ORDER,
            ):
                handle.write(line + "\n")
            handle.write("\n")

        if has_derived and biology_counts:
            section += 1
            handle.write(f"{section}. SHARED-LOCUS BIOLOGY (BOTH CLADES PRESENT)\n")
            handle.write("-" * 40 + "\n")
            biology_order = [
                "shared_ancestral_fixed_both_clades",
                "shared_ancestral_not_fixed_both_clades",
                "shared_independent_insertion",
                "shared_ambiguous_cross_group",
                "shared_other_cross_group",
                "within_clade_only",
            ]
            biology_labels = {
                "shared_ancestral_fixed_both_clades": (
                    "ancestral + fixed_present in both clades"
                ),
                "shared_ancestral_not_fixed_both_clades": (
                    "ancestral but not fixed in both clades"
                ),
                "shared_independent_insertion": "independent cross-group insertion",
                "shared_ambiguous_cross_group": "ambiguous cross-group origin",
                "shared_other_cross_group": "other cross-group origin",
                "within_clade_only": "present in only one clade",
            }
            for key in biology_order:
                count = biology_counts.get(key, 0)
                if count == 0 and key.startswith("shared_other"):
                    continue
                pct = 100 * count / total_groups if total_groups else 0
                handle.write(
                    f"  {biology_labels[key]:<42s}  {count:>7,d}  ({pct:5.1f}%)\n"
                )
            handle.write("\n")

        section += 1
        handle.write(f"{section}. GENE OCCUPANCY\n")
        handle.write("-" * 40 + "\n")
        handle.write("Note: gene annotations from validation step; may be unreliable.\n")
        if g1_gene_total:
            handle.write(
                f"Group 1: {g1_in_gene:,}/{g1_gene_total:,} "
                f"({100 * g1_in_gene / g1_gene_total:.1f}%) within annotated genes\n"
            )
        if g2_gene_total:
            handle.write(
                f"Group 2: {g2_in_gene:,}/{g2_gene_total:,} "
                f"({100 * g2_in_gene / g2_gene_total:.1f}%) within annotated genes\n"
            )
        handle.write("\n")

        section += 1
        handle.write(f"{section}. FAMILY FREQUENCY COMPOSITION\n")
        handle.write("-" * 40 + "\n")
        all_families = sorted(set(g1_fam_counts.keys()) | set(g2_fam_counts.keys()))
        g1_fam_total = sum(g1_fam_counts.values())
        g2_fam_total = sum(g2_fam_counts.values())
        handle.write(
            f"{'Family':>8s}  {'Group1':>8s} {'(%)':>7s}  {'Group2':>8s} {'(%)':>7s}\n"
        )
        for fam in all_families:
            g1_c = g1_fam_counts.get(fam, 0)
            g2_c = g2_fam_counts.get(fam, 0)
            g1_pct = 100 * g1_c / g1_fam_total if g1_fam_total else 0
            g2_pct = 100 * g2_c / g2_fam_total if g2_fam_total else 0
            handle.write(
                f"{int(fam):>8d}  {g1_c:>8,} {g1_pct:>6.1f}%  "
                f"{g2_c:>8,} {g2_pct:>6.1f}%\n"
            )
        handle.write(
            f"{'Total':>8s}  {g1_fam_total:>8,}          {g2_fam_total:>8,}\n"
        )
        handle.write("\n" + "=" * 72 + "\n")

    print(f"Summary written to: {args.output}")


if __name__ == "__main__":
    main()
