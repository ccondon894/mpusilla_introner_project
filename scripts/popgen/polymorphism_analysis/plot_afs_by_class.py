"""Comparative derived AFS for synonymous SNPs, nonsynonymous SNPs, and
polymorphic introners. SNPs are polarized using G2 as the outgroup.

Inputs:
    --snpeff_vcf       SnpEff-annotated, MT-removed VCF (text)
    --genotype_matrix  genotype_matrix.final.tsv (introner data)
    --non_introner_matrix  intron_genotype_matrix.tsv (non-introner intron data)
    --group1           Comma-separated G1 sample list (n=11 typical)
    --outgroup         Comma-separated G2 sample list (used to polarize SNPs)
    --output_tsv       Per-class SFS counts (long-format TSV)
    --output_plot      Comparative density bar plot

Outputs a comparative G1 AFS (length n-1, count bins 1..n-1) for four classes:
introner present counts, non-introner intron present counts, synonymous derived
SNP counts, and nonsynonymous derived SNP counts. Each class is normalized to a
density (proportion within that class).

Genotype matrix presence convention:
    presence == 1  → present
    presence == 2  → absent
    presence == 3  → missing
"""

import argparse
import gzip
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

SCRIPTS_DIR = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(SCRIPTS_DIR))

from figure_color_guide import DEFAULT_GUIDE_PATH, get_color, load_color_guide


SYN_TERMS = {"synonymous_variant"}
NONSYN_TERMS = {"missense_variant"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--snpeff_vcf", required=True)
    p.add_argument("--genotype_matrix", required=True)
    p.add_argument("--non_introner_matrix", required=True)
    p.add_argument("--group1", required=True, help="Comma-separated G1 sample list")
    p.add_argument("--outgroup", required=True, help="Comma-separated G2 sample list")
    p.add_argument("--output_tsv", required=True)
    p.add_argument("--output_pdf", required=True)
    p.add_argument("--output_png", required=True)
    p.add_argument("--color_guide", default=str(DEFAULT_GUIDE_PATH))
    p.add_argument(
        "--status_filter",
        default="consistent,singleton",
        help="within_group_status values to keep for introners",
    )
    return p.parse_args()


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_consequence(ann_field):
    """Return the first ANN entry's consequence (col 1 after the allele)."""
    if not ann_field:
        return None
    first = ann_field.split(",")[0]
    parts = first.split("|")
    return parts[1] if len(parts) > 1 else None


def haploid_call(gt_str):
    """Return 0/1 for haploid genotype, or None if missing/heterozygous."""
    if not gt_str or gt_str == "." or gt_str.startswith("."):
        return None
    g = gt_str.split(":")[0].split("/")[0].split("|")[0]
    if g == ".":
        return None
    return int(g) if g in ("0", "1") else None


def snp_sfs_by_class(vcf_path, g1_idx, g2_idx, n_g1):
    """Build per-class SFS arrays of length n_g1+1 (folded later to 1..n-1).

    Returns: dict { 'synonymous': counts, 'nonsynonymous': counts, ... }.
    Polarization rule: derived = allele not seen in G2 (G2 must be unanimous
    and fully called); sites with split or missing G2 are skipped.
    """
    syn = np.zeros(n_g1 + 1, dtype=np.int64)
    nonsyn = np.zeros(n_g1 + 1, dtype=np.int64)
    drops = Counter()

    with open_text(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue

            # Bi-allelic SNPs only.
            if len(fields[3]) != 1 or len(fields[4]) != 1 or "," in fields[4]:
                drops["non_biallelic"] += 1
                continue

            info = fields[7]
            ann = None
            for kv in info.split(";"):
                if kv.startswith("ANN="):
                    ann = kv[4:]
                    break
            consequence = parse_consequence(ann)
            if consequence in SYN_TERMS:
                target = syn
            elif consequence in NONSYN_TERMS:
                target = nonsyn
            else:
                drops["non_coding_or_other"] += 1
                continue

            samples = fields[9:]

            g2_calls = [haploid_call(samples[i]) for i in g2_idx]
            if any(c is None for c in g2_calls):
                drops["g2_missing"] += 1
                continue
            g2_set = set(g2_calls)
            if len(g2_set) != 1:
                drops["g2_polymorphic"] += 1
                continue
            ancestral = g2_set.pop()
            derived = 1 - ancestral

            g1_calls = [haploid_call(samples[i]) for i in g1_idx]
            if any(c is None for c in g1_calls):
                drops["g1_missing"] += 1
                continue
            derived_count = sum(1 for c in g1_calls if c == derived)
            target[derived_count] += 1

    print(
        f"SNP drops: {dict(drops)}",
        file=sys.stderr,
    )
    return {"synonymous": syn, "nonsynonymous": nonsyn}


def vcf_sample_indices(vcf_path, g1_samples, g2_samples):
    """Return (g1_idx, g2_idx) with positions in the VCF sample columns."""
    with open_text(vcf_path) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")[9:]
                idx = {s: i for i, s in enumerate(cols)}
                missing = [s for s in g1_samples + g2_samples if s not in idx]
                if missing:
                    raise SystemExit(f"Samples missing from VCF: {missing}")
                return [idx[s] for s in g1_samples], [idx[s] for s in g2_samples]
    raise SystemExit("No #CHROM header line found")


def introner_sfs(matrix_path, g1_samples, status_filter):
    """Build the introner present-allele-count SFS for G1."""
    n = len(g1_samples)
    keep_status = set(s.strip() for s in status_filter.split(",") if s.strip())

    df = pd.read_csv(matrix_path, sep="\t")
    df = df[df["sample"].isin(g1_samples)].copy()
    pivot = df.pivot(index="ortholog_id", columns="sample", values="presence")
    missing_cols = [s for s in g1_samples if s not in pivot.columns]
    if missing_cols:
        raise SystemExit(f"Samples missing from genotype matrix: {missing_cols}")
    pivot = pivot[g1_samples]

    status = (
        df.groupby("ortholog_id")["within_group_status"].first().reindex(pivot.index)
    )
    pivot = pivot[status.isin(keep_status)]

    has_missing = (pivot == 3).any(axis=1)
    pivot = pivot[~has_missing]

    present_counts = (pivot == 1).sum(axis=1).astype(int)
    sfs_full = np.zeros(n + 1, dtype=np.int64)
    for k in range(n + 1):
        sfs_full[k] = int((present_counts == k).sum())
    return sfs_full


def non_introner_intron_sfs(matrix_path, g1_samples):
    """Build the non-introner intron present-allele-count SFS for G1."""
    n = len(g1_samples)

    df = pd.read_csv(matrix_path, sep="\t")
    missing_cols = [s for s in g1_samples if s not in df.columns]
    if missing_cols:
        raise SystemExit(f"Samples missing from non-introner matrix: {missing_cols}")

    sub = df[g1_samples].copy()
    has_missing = (sub == 3).any(axis=1)
    sub = sub[~has_missing]

    present_counts = (sub == 1).sum(axis=1).astype(int)
    sfs_full = np.zeros(n + 1, dtype=np.int64)
    for k in range(n + 1):
        sfs_full[k] = int((present_counts == k).sum())
    return sfs_full


def main():
    args = parse_args()
    color_guide = load_color_guide(args.color_guide)
    class_colors = {
        "synonymous": get_color("Synonymous mutation", color_guide),
        "nonsynonymous": get_color("Nonsynonymous mutation", color_guide),
        "introner": get_color("Introner", color_guide),
        "non_introner_intron": get_color("Intron", color_guide),
    }

    g1 = [s.strip() for s in args.group1.split(",") if s.strip()]
    g2 = [s.strip() for s in args.outgroup.split(",") if s.strip()]
    n = len(g1)

    print(f"G1 (n={n}): {g1}", file=sys.stderr)
    print(f"G2 (n={len(g2)}): {g2}", file=sys.stderr)

    g1_idx, g2_idx = vcf_sample_indices(args.snpeff_vcf, g1, g2)
    snp_sfs = snp_sfs_by_class(args.snpeff_vcf, g1_idx, g2_idx, n)

    intr_sfs = introner_sfs(args.genotype_matrix, g1, args.status_filter)
    non_intr_sfs = non_introner_intron_sfs(args.non_introner_matrix, g1)

    seg = slice(1, n)  # bins 1..n-1, exclude fixed
    intr_seg = intr_sfs[seg].astype(float)
    non_intr_seg = non_intr_sfs[seg].astype(float)
    syn_seg = snp_sfs["synonymous"][seg].astype(float)
    nonsyn_seg = snp_sfs["nonsynonymous"][seg].astype(float)

    def density(arr):
        total = arr.sum()
        return arr / total if total > 0 else arr

    intr_d = density(intr_seg)
    non_intr_d = density(non_intr_seg)
    syn_d = density(syn_seg)
    nonsyn_d = density(nonsyn_seg)

    bins = np.arange(1, n)
    rows = []
    for i, b in enumerate(bins):
        rows.append({
            "derived_count": int(b),
            "introner_count": int(intr_seg[i]),
            "non_introner_intron_count": int(non_intr_seg[i]),
            "synonymous_count": int(syn_seg[i]),
            "nonsynonymous_count": int(nonsyn_seg[i]),
            "introner_density": float(intr_d[i]),
            "non_introner_intron_density": float(non_intr_d[i]),
            "synonymous_density": float(syn_d[i]),
            "nonsynonymous_density": float(nonsyn_d[i]),
        })
    pd.DataFrame(rows).to_csv(args.output_tsv, sep="\t", index=False)

    width = 0.20
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(bins - 1.5 * width, syn_d, width,
           label="Synonymous",
           color=class_colors["synonymous"], edgecolor="black")
    ax.bar(bins - 0.5 * width, nonsyn_d, width,
           label="Nonsynonymous",
           color=class_colors["nonsynonymous"], edgecolor="black")
    ax.bar(bins + 0.5 * width, intr_d, width,
           label="Introners",
           color=class_colors["introner"], edgecolor="black")
    ax.bar(bins + 1.5 * width, non_intr_d, width,
           label="Introns",
           color=class_colors["non_introner_intron"], edgecolor="black")
    ax.set_xticks(bins)
    ax.set_xlabel("Group 1 allele count (derived SNPs / present introns)")
    ax.set_ylabel("Density")
    ax.legend(loc="upper right", bbox_to_anchor=(0.985, 0.985), frameon=False)
    fig.tight_layout()
    fig.savefig(args.output_pdf, dpi=300, bbox_inches="tight", format="pdf")
    fig.savefig(args.output_png, dpi=300, bbox_inches="tight", format="png")
    plt.close(fig)

    print(
        f"Introner total seg: {int(intr_seg.sum())}  "
        f"Non-introner intron total seg: {int(non_intr_seg.sum())}  "
        f"Syn total: {int(syn_seg.sum())}  Nonsyn total: {int(nonsyn_seg.sum())}",
        file=sys.stderr,
    )
    print(f"Wrote {args.output_tsv}, {args.output_pdf}", file=sys.stderr)


if __name__ == "__main__":
    main()
