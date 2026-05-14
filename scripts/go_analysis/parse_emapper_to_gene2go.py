#!/usr/bin/env python3
"""
Parse emapper.py annotations output into a gene2go.json compatible with the
existing GO enrichment pipeline.

Protein FASTA headers look like:
    >182920.3.0.228 gene=gm1.1_g.3.0.228 seq_id=scaffold_1 type=cds

emapper uses the first whitespace-delimited token as the query name (e.g.
"182920.3.0.228"), which is the transcriptName + ".3.0.228" suffix. The
annotation file's transcriptName and peptideName diverge for some genes,
so we re-parse the FASTA headers to build a query -> locusName map
directly from the `gene=` field.
"""

import argparse
import json
import re
from pathlib import Path

import pandas as pd


def parse_emapper_annotations(emapper_tsv: Path) -> pd.DataFrame:
    """Read emapper .annotations file with its comment-prefixed header."""
    with open(emapper_tsv) as handle:
        header_line = None
        for line in handle:
            if line.startswith("#query"):
                header_line = line.lstrip("#").rstrip("\n")
                break
        if header_line is None:
            raise RuntimeError(f"No '#query' header found in {emapper_tsv}")
        cols = header_line.split("\t")
    return pd.read_csv(emapper_tsv, sep="\t", comment="#", names=cols, dtype=str)


def build_query_to_locus(fasta_path: Path) -> dict[str, str]:
    """Build query_id -> locusName map from protein FASTA headers."""
    q2l: dict[str, str] = {}
    suffix_re = re.compile(r"\.\d+\.\d+\.\d+$")
    with open(fasta_path) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            line = line[1:].rstrip("\n")
            tokens = line.split()
            query = tokens[0]
            locus = None
            for token in tokens[1:]:
                if token.startswith("gene="):
                    locus = suffix_re.sub("", token[5:])
                    break
            if locus is not None:
                q2l[query] = locus
    return q2l


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--emapper", required=True, help="emapper .annotations TSV")
    parser.add_argument("--fasta", required=True, help="Reference proteins FASTA used as emapper input")
    parser.add_argument("--annotation", required=True, help="Annotation info TSV used for the gene universe")
    parser.add_argument("--gene2go", required=True, help="Output gene2go.json")
    parser.add_argument("--stats", required=True, help="Output coverage stats")
    args = parser.parse_args()

    emapper_df = parse_emapper_annotations(Path(args.emapper))
    query_to_locus = build_query_to_locus(Path(args.fasta))
    emapper_df["locusName"] = emapper_df["query"].map(query_to_locus)
    missing = emapper_df["locusName"].isna().sum()

    annotation_df = pd.read_csv(args.annotation, sep="\t", dtype=str, na_filter=False)

    gene2go: dict[str, set[str]] = {}
    for _, row in emapper_df.iterrows():
        locus = row["locusName"]
        if pd.isna(locus):
            continue
        go_field = row.get("GOs", "")
        if not go_field or go_field == "-":
            continue
        terms = {term for term in go_field.split(",") if term.startswith("GO:")}
        if not terms:
            continue
        gene2go.setdefault(locus, set()).update(terms)

    all_loci = set(annotation_df["locusName"])
    output = {locus: sorted(gene2go.get(locus, set())) for locus in all_loci}
    Path(args.gene2go).write_text(json.dumps(output, indent=2))

    n_genes = len(output)
    n_with_go = sum(1 for values in output.values() if values)
    total_go = sum(len(values) for values in output.values())
    with open(args.stats, "w") as handle:
        handle.write("eggNOG-mapper GO Coverage Statistics\n")
        handle.write("=" * 40 + "\n")
        handle.write(f"Total genes (annotation universe): {n_genes}\n")
        handle.write(f"emapper queries (proteins hit): {len(emapper_df)}\n")
        handle.write(f"emapper queries with no locus match: {missing}\n")
        handle.write(f"Genes with GO terms: {n_with_go} ({100 * n_with_go / n_genes:.1f}%)\n")
        handle.write(f"Total GO annotations: {total_go}\n")
        if n_with_go:
            handle.write(f"Average GO terms per annotated gene: {total_go / n_with_go:.1f}\n")

    print(f"Wrote {n_genes} genes ({n_with_go} with GO) -> {args.gene2go}")


if __name__ == "__main__":
    main()
