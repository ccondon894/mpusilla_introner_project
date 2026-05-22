#!/usr/bin/env python3
"""
Summarize subgroup-specific GO terms relative to the all-introner signal.

The input is the empirically annotated introner-group GO table. For each
empirically significant subgroup GO row, this script reports whether the GO ID
is absent from the all_introners significant set and whether it is also absent
after accounting for GO parent terms already present in all_introners.
"""

import argparse
import sys
from functools import lru_cache

import pandas as pd

from introner_group_go_enrichment import parse_obo


def parse_args():
    parser = argparse.ArgumentParser(
        description="Identify introner subgroup GO terms not captured by all_introners."
    )
    parser.add_argument("--empirical-go", required=True)
    parser.add_argument("--go-obo", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--unique-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--base-group", default="all_introners")
    parser.add_argument("--significance-column", default="empirical_significant_global")
    parser.add_argument("--fdr-column", default="empirical_fdr_bh_global")
    return parser.parse_args()


def format_ancestor_hits(hits, terms):
    labels = []
    for go_id in sorted(hits):
        name = terms.get(go_id, {}).get("name", "")
        labels.append(f"{go_id} ({name})" if name else go_id)
    return ";".join(labels)


def write_summary(path, annotated, unique_rows, args):
    with open(path, "w") as handle:
        handle.write("Introner Group Unique GO Term Summary\n")
        handle.write("=" * 45 + "\n\n")
        handle.write(f"Empirical GO table: {args.empirical_go}\n")
        handle.write(f"Base group: {args.base_group}\n")
        handle.write(f"Significance column: {args.significance_column}\n")
        handle.write(
            "Exact unique: subgroup GO ID is absent from significant base-group GO IDs.\n"
        )
        handle.write(
            "Ontology unique: exact unique and no subgroup GO ancestor is significant in the base group.\n\n"
        )

        base_count = int((annotated["introner_group"] == args.base_group).sum())
        subgroup = annotated[annotated["introner_group"] != args.base_group]
        handle.write(f"Significant base-group GO rows: {base_count}\n")
        handle.write(f"Significant subgroup GO rows evaluated: {len(subgroup)}\n")
        handle.write(f"Ontology-unique subgroup GO rows: {len(unique_rows)}\n\n")

        for group in sorted(subgroup["introner_group"].unique()):
            rows = subgroup[subgroup["introner_group"] == group]
            exact = int(rows["exact_unique_vs_base"].sum())
            ontology = int(rows["ontology_unique_vs_base"].sum())
            handle.write(
                f"{group}: {len(rows)} significant rows, "
                f"{exact} exact unique, {ontology} ontology unique\n"
            )
            group_unique = unique_rows[unique_rows["introner_group"] == group]
            for _, row in group_unique.sort_values(
                ["enrichment_type", args.fdr_column, "GO_term"]
            ).iterrows():
                handle.write(
                    f"  {row['GO_term']}\t{row['enrichment_type']}\t"
                    f"{row['term_name']}\t{args.fdr_column}={row[args.fdr_column]:.6g}\n"
                )
            handle.write("\n")


def main():
    args = parse_args()

    print("Loading empirical GO table", file=sys.stderr)
    df = pd.read_csv(args.empirical_go, sep="\t")
    required = {
        "introner_group",
        "GO_term",
        "term_name",
        "enrichment_type",
        args.significance_column,
        args.fdr_column,
    }
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError("Input table missing columns: " + ", ".join(missing))

    print("Loading GO ontology", file=sys.stderr)
    terms = parse_obo(args.go_obo)

    significant = df[df[args.significance_column].astype(bool)].copy()
    if args.base_group not in set(significant["introner_group"]):
        raise ValueError(
            f"Base group {args.base_group!r} has no significant GO rows in input table."
        )

    base_terms = set(
        significant.loc[significant["introner_group"] == args.base_group, "GO_term"]
    )

    @lru_cache(maxsize=None)
    def ancestors(go_id):
        seen = set()
        stack = list(terms.get(go_id, {}).get("parents", []))
        while stack:
            parent = stack.pop()
            if parent in seen:
                continue
            seen.add(parent)
            stack.extend(terms.get(parent, {}).get("parents", []))
        return frozenset(seen)

    significant["base_group"] = args.base_group
    significant["exact_unique_vs_base"] = ~significant["GO_term"].isin(base_terms)
    significant.loc[
        significant["introner_group"] == args.base_group,
        "exact_unique_vs_base",
    ] = False
    significant["base_ancestor_hits"] = significant["GO_term"].map(
        lambda go_id: format_ancestor_hits(ancestors(go_id) & base_terms, terms)
    )
    significant["ontology_unique_vs_base"] = (
        significant["exact_unique_vs_base"] & significant["base_ancestor_hits"].eq("")
    )

    significant = significant.sort_values(
        [
            "introner_group",
            "ontology_unique_vs_base",
            "exact_unique_vs_base",
            args.fdr_column,
            "GO_term",
        ],
        ascending=[True, False, False, True, True],
    )
    unique_rows = significant[
        (significant["introner_group"] != args.base_group)
        & significant["ontology_unique_vs_base"]
    ].copy()

    significant.to_csv(args.output, sep="\t", index=False, float_format="%.6g")
    unique_rows.to_csv(args.unique_output, sep="\t", index=False, float_format="%.6g")
    write_summary(args.summary_output, significant, unique_rows, args)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
