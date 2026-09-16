#!/usr/bin/env python3
"""Rank genes with clean introners, CCMP1545/RCC1614 polymorphism, and R2C2 isoforms.

Gene mapping is rebuilt from the current genotype matrix and CCMP1545 GTF each
run because ortholog_id values can be reassigned when the genotype matrix is
regenerated.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


PROJECT = Path("/scratch1/chris/mpusilla_introner_project")
SENSITIVE = Path(
    "/scratch1/chris/introner-expression-analysis/R2C2_data/sensitiveIsoforms"
)
GENOTYPE = PROJECT / "results/genotyping/genotype_matrix.final.tsv"
CCMP1545_ANNOTATION_GTF = PROJECT / "results/annotations/CCMP1545.gtf"
CURRENT_INTRONER_INTRON_MAP = (
    PROJECT
    / "analysis/ccmp1545_gtf_liftover_test/current_genotype_introner_vs_annotated_intron_size.tsv"
)
CCMP1545_GTF = SENSITIVE / "09092025_834_Isoforms.filtered.clean.scaffold_corrected.sorted.gtf"
RCC1614_GTF = SENSITIVE / "09092025_1614_Isoforms.filtered.clean.sorted.gtf"
CCMP1545_SQANTI = SENSITIVE / "sqanti3_output_834/834_isoforms_classification.filtered.txt"
RCC1614_SQANTI = SENSITIVE / "sqanti3_output_1614/1614_isoforms_classification.filtered.txt"
OUT = (
    PROJECT
    / "analysis/ccmp1545_gtf_liftover_test/pairwise_1614_isoform_gene_candidates.tsv"
)
GENOTYPE_FLANK = 100
MIN_BODY_IN_INTRON_FRAC = 0.95


def parse_attrs(attr_text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in attr_text.strip().rstrip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if " " in item:
            key, value = item.split(" ", 1)
            attrs[key] = value.strip().strip('"')
        elif "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value.strip().strip('"')
    return attrs


def norm_gene(gene: object) -> str:
    value = "" if pd.isna(gene) else str(gene)
    if value.endswith(".3.0.228"):
        value = value[: -len(".3.0.228")]
    if value.startswith("novelGene_") and value.endswith("_AS"):
        value = value[len("novelGene_") : -len("_AS")]
        if value.endswith(".3.0.228"):
            value = value[: -len(".3.0.228")]
    return value


def transcript_counts(gtf_path: Path) -> tuple[Counter[str], Counter[str]]:
    transcripts_by_gene: dict[str, set[str]] = defaultdict(set)
    exons_by_transcript: dict[tuple[str, str], int] = defaultdict(int)
    with gtf_path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            feature = fields[2]
            attrs = parse_attrs(fields[8])
            gene = norm_gene(attrs.get("gene_id") or attrs.get("gene_name"))
            tx = attrs.get("transcript_id")
            if not gene or not tx:
                continue
            if feature == "transcript":
                transcripts_by_gene[gene].add(tx)
            elif feature == "exon":
                exons_by_transcript[(gene, tx)] += 1

    isoform_count = Counter({gene: len(txs) for gene, txs in transcripts_by_gene.items()})
    multiexon_count: Counter[str] = Counter()
    for gene, txs in transcripts_by_gene.items():
        multiexon_count[gene] = sum(1 for tx in txs if exons_by_transcript[(gene, tx)] >= 2)
    return isoform_count, multiexon_count


def sqanti_summary(path: Path) -> tuple[Counter[str], Counter[str], dict[str, Counter[str]], dict[str, Counter[str]]]:
    if not path.exists():
        return Counter(), Counter(), {}, {}
    df = pd.read_csv(path, sep="\t")
    gene_col = "associated_gene" if "associated_gene" in df.columns else "Gene"
    df["norm_gene"] = df[gene_col].map(norm_gene)
    retained = df["subcategory"].fillna("").str.lower().str.contains("intron_retention")
    retained_count = Counter(df.loc[retained, "norm_gene"])
    isoform_count = Counter(df["norm_gene"])
    categories = {
        gene: Counter(sub["structural_category"].fillna(""))
        for gene, sub in df.groupby("norm_gene")
    }
    subcategories = {
        gene: Counter(sub["subcategory"].fillna(""))
        for gene, sub in df.groupby("norm_gene")
    }
    return isoform_count, retained_count, categories, subcategories


def compact_counts(counter: Counter[str], limit: int = 4) -> str:
    return ";".join(f"{key}:{value}" for key, value in counter.most_common(limit) if key)


def read_reference_introns(gtf_path: Path) -> dict[str, list[dict[str, object]]]:
    exons_by_tx: dict[tuple[str, str], list[tuple[str, int, int]]] = defaultdict(list)
    tx_gene: dict[tuple[str, str], str] = {}
    with gtf_path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            attrs = parse_attrs(fields[8])
            gene = attrs.get("gene_id") or attrs.get("gene_name")
            tx = attrs.get("transcript_id")
            if not gene or not tx:
                continue
            key = (fields[0], tx)
            tx_gene[key] = gene
            exons_by_tx[key].append((fields[0], int(fields[3]), int(fields[4])))

    introns_by_contig: dict[str, list[dict[str, object]]] = defaultdict(list)
    for key, exons in exons_by_tx.items():
        contig, tx = key
        gene = tx_gene[key]
        exons_sorted = sorted(exons, key=lambda item: (item[1], item[2]))
        for intron_index, (left, right) in enumerate(
            zip(exons_sorted, exons_sorted[1:]), start=1
        ):
            intron_start = left[2] + 1
            intron_end = right[1] - 1
            if intron_start > intron_end:
                continue
            introns_by_contig[contig].append(
                {
                    "gene": gene,
                    "transcript": tx,
                    "intron_index": intron_index,
                    "intron_start": intron_start,
                    "intron_end": intron_end,
                    "intron_len": intron_end - intron_start + 1,
                }
            )
    for contig in introns_by_contig:
        introns_by_contig[contig].sort(key=lambda item: (item["intron_start"], item["intron_end"]))
    return introns_by_contig


def interval_overlap(
    start_a: int,
    end_a: int,
    start_b: int,
    end_b: int,
) -> int:
    start = max(start_a, start_b)
    end = min(end_a, end_b)
    return max(0, end - start + 1)


def build_current_clean_introner_map(
    genotype: pd.DataFrame,
    gtf_path: Path,
    flank: int = GENOTYPE_FLANK,
    min_body_in_intron_frac: float = MIN_BODY_IN_INTRON_FRAC,
) -> pd.DataFrame:
    introns_by_contig = read_reference_introns(gtf_path)
    ccmp_present = genotype[
        genotype["sample"].eq("CCMP1545")
        & genotype["presence"].eq(1)
        & genotype["start"].ge(0)
        & genotype["end"].ge(0)
    ].copy()

    rows = []
    for _, locus in ccmp_present.iterrows():
        contig = str(locus["contig"])
        body_start = int(locus["start"]) + flank
        body_end = int(locus["end"]) - flank
        if body_end < body_start:
            continue
        body_len = body_end - body_start + 1
        best = None
        for intron in introns_by_contig.get(contig, []):
            overlap = interval_overlap(
                body_start,
                body_end,
                int(intron["intron_start"]),
                int(intron["intron_end"]),
            )
            if overlap == 0:
                continue
            body_frac = overlap / body_len
            intron_frac = overlap / int(intron["intron_len"])
            length_ratio = body_len / int(intron["intron_len"])
            candidate = {
                "ortholog_id": locus["ortholog_id"],
                "gene": intron["gene"],
                "contig": contig,
                "body_start": body_start,
                "body_end": body_end,
                "body_len": body_len,
                "transcript": intron["transcript"],
                "intron_index": intron["intron_index"],
                "intron_start": intron["intron_start"],
                "intron_end": intron["intron_end"],
                "intron_len": intron["intron_len"],
                "overlap_bp": overlap,
                "body_in_intron_frac": body_frac,
                "intron_covered_frac": intron_frac,
                "length_ratio": length_ratio,
                "group1_pattern": locus["group1_pattern"],
                "family": locus["family"],
            }
            score = (body_frac, overlap, intron_frac)
            if best is None or score > best[0]:
                best = (score, candidate)
        if best is not None and best[1]["body_in_intron_frac"] >= min_body_in_intron_frac:
            rows.append(best[1])

    if not rows:
        return pd.DataFrame(
            columns=[
                "ortholog_id",
                "gene",
                "contig",
                "body_start",
                "body_end",
                "body_len",
                "transcript",
                "intron_index",
                "intron_start",
                "intron_end",
                "intron_len",
                "overlap_bp",
                "body_in_intron_frac",
                "intron_covered_frac",
                "length_ratio",
                "group1_pattern",
                "family",
            ]
        )
    clean = pd.DataFrame(rows)
    clean["norm_gene"] = clean["gene"].map(norm_gene)
    return clean.sort_values(["norm_gene", "ortholog_id"]).reset_index(drop=True)


def main() -> None:
    global PROJECT, SENSITIVE, GENOTYPE, CCMP1545_ANNOTATION_GTF, CURRENT_INTRONER_INTRON_MAP
    global CCMP1545_GTF, RCC1614_GTF, CCMP1545_SQANTI, RCC1614_SQANTI, OUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--isoform-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sidecar", type=Path, required=True)
    args = parser.parse_args()
    PROJECT, SENSITIVE, OUT = args.project_root, args.isoform_dir, args.output
    GENOTYPE = PROJECT / "results/genotyping/genotype_matrix.final.tsv"
    CCMP1545_ANNOTATION_GTF = PROJECT / "results/annotations/CCMP1545.gtf"
    CURRENT_INTRONER_INTRON_MAP = args.sidecar
    CCMP1545_GTF = SENSITIVE / "09092025_834_Isoforms.filtered.clean.scaffold_corrected.sorted.gtf"
    RCC1614_GTF = SENSITIVE / "09092025_1614_Isoforms.filtered.clean.sorted.gtf"
    CCMP1545_SQANTI = PROJECT / "data/sqanti3_output_834/834_isoforms_classification.filtered.txt"
    RCC1614_SQANTI = PROJECT / "data/sqanti3_output_1614/1614_isoforms_classification.filtered.txt"
    genotype = pd.read_csv(GENOTYPE, sep="\t")

    clean = build_current_clean_introner_map(genotype, CCMP1545_ANNOTATION_GTF)
    CURRENT_INTRONER_INTRON_MAP.parent.mkdir(parents=True, exist_ok=True)
    clean.to_csv(CURRENT_INTRONER_INTRON_MAP, sep="\t", index=False)

    ccmp_iso, ccmp_multiexon = transcript_counts(CCMP1545_GTF)
    rcc_iso, rcc_multiexon = transcript_counts(RCC1614_GTF)
    _, ccmp_retained, ccmp_categories, ccmp_subcategories = sqanti_summary(CCMP1545_SQANTI)
    _, rcc_retained, rcc_categories, rcc_subcategories = sqanti_summary(RCC1614_SQANTI)

    pair = genotype[
        genotype["sample"].isin(["CCMP1545", "RCC1614"])
        & genotype["presence"].isin([1, 2])
    ].copy()
    pair_matrix = pair.pivot_table(
        index="ortholog_id", columns="sample", values="presence", aggfunc="first"
    )
    pair_polymorphic = pair_matrix[
        pair_matrix.get("CCMP1545").isin([1, 2])
        & pair_matrix.get("RCC1614").isin([1, 2])
        & (pair_matrix["CCMP1545"] != pair_matrix["RCC1614"])
    ]
    pair_oids = set(pair_polymorphic.index)

    rows = []
    for gene, sub in clean.groupby("norm_gene"):
        clean_oids = sorted(set(sub["ortholog_id"]))
        pair_clean_oids = [oid for oid in clean_oids if oid in pair_oids]
        if len(clean_oids) < 2 or not pair_clean_oids:
            continue
        if ccmp_iso[gene] < 2 or rcc_iso[gene] < 2:
            continue

        pair_patterns = []
        for oid in pair_clean_oids:
            calls = pair_matrix.loc[oid]
            pair_patterns.append(
                f"{oid}:CCMP1545={int(calls['CCMP1545'])},RCC1614={int(calls['RCC1614'])}"
            )

        group1_patterns = Counter(sub["group1_pattern"])
        score = (
            len(clean_oids) * 10
            + len(pair_clean_oids) * 20
            + min(ccmp_iso[gene], 8)
            + min(rcc_iso[gene], 8)
            + 3 * min(ccmp_retained[gene], 3)
            + 3 * min(rcc_retained[gene], 3)
            + 5 * float(sub["body_in_intron_frac"].median())
        )
        rows.append(
            {
                "gene": gene,
                "score": score,
                "clean_introner_count": len(clean_oids),
                "ccmp1545_rcc1614_polymorphic_clean_count": len(pair_clean_oids),
                "ccmp1545_isoform_count": ccmp_iso[gene],
                "rcc1614_isoform_count": rcc_iso[gene],
                "ccmp1545_multiexon_isoform_count": ccmp_multiexon[gene],
                "rcc1614_multiexon_isoform_count": rcc_multiexon[gene],
                "ccmp1545_retained_intron_isoform_count": ccmp_retained[gene],
                "rcc1614_retained_intron_isoform_count": rcc_retained[gene],
                "median_body_host_ratio": sub["length_ratio"].median(),
                "median_body_in_intron_frac": sub["body_in_intron_frac"].median(),
                "median_host_intron_covered_frac": sub["intron_covered_frac"].median(),
                "group1_patterns": compact_counts(group1_patterns),
                "pairwise_polymorphic_clean_orthologs": ",".join(pair_clean_oids),
                "pairwise_calls": ";".join(pair_patterns),
                "clean_orthologs": ",".join(clean_oids),
                "ccmp1545_sqanti_categories": compact_counts(ccmp_categories.get(gene, Counter())),
                "rcc1614_sqanti_categories": compact_counts(rcc_categories.get(gene, Counter())),
                "ccmp1545_sqanti_subcategories": compact_counts(ccmp_subcategories.get(gene, Counter())),
                "rcc1614_sqanti_subcategories": compact_counts(rcc_subcategories.get(gene, Counter())),
            }
        )

    rows.sort(
        key=lambda row: (
            row["score"],
            row["ccmp1545_rcc1614_polymorphic_clean_count"],
            row["clean_introner_count"],
            min(row["ccmp1545_isoform_count"], row["rcc1614_isoform_count"]),
        ),
        reverse=True,
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]) if rows else ["gene"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {OUT}")
    print(f"Candidates: {len(rows)}")
    for row in rows[:10]:
        print(
            "\t".join(
                [
                    row["gene"],
                    f"clean={row['clean_introner_count']}",
                    f"pair_poly={row['ccmp1545_rcc1614_polymorphic_clean_count']}",
                    f"isoforms={row['ccmp1545_isoform_count']}/{row['rcc1614_isoform_count']}",
                    f"retained={row['ccmp1545_retained_intron_isoform_count']}/{row['rcc1614_retained_intron_isoform_count']}",
                    row["pairwise_polymorphic_clean_orthologs"],
                ]
            )
        )


if __name__ == "__main__":
    main()
