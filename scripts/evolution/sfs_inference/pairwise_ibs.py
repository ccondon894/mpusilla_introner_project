"""Pairwise IBS / Hamming-distance diagnostic for cluster_1.

Tests the hypothesis that 2 of the 6 cluster_1 isolates are near-clonal,
which would produce a bin-2 bump in the folded SFS via "fixed in the pair,
absent in the other four" sites. Outputs:

  - pairwise IBS matrix (fraction of segregating sites with identical alleles)
  - heatmap PNG
  - per-pair private-allele decomposition of the bin-2 (AC=2 or AC=4) sites:
    is the bin-2 mass dominated by one specific pair, or spread evenly?
"""

import argparse
import gzip
import sys
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf", required=True)
    p.add_argument("--samples", required=True, help="Comma-separated sample list")
    p.add_argument("--ibs_matrix", required=True, help="TSV — full IBS matrix")
    p.add_argument("--ibs_summary", required=True, help="TSV — pairwise IBS list sorted by IBS")
    p.add_argument("--bin2_decomp", required=True, help="TSV — bin-2 ownership counts")
    p.add_argument("--plot", required=True, help="PNG with heatmap + bar of pair-fixed counts")
    return p.parse_args()


def _open(path):
    return gzip.open(path, "rt") if path.endswith((".gz", ".bgz")) else open(path)


def build_genotype_matrix(vcf_path, sample_list):
    n = len(sample_list)
    sample_cols = None
    rows = []
    n_skipped_missing = 0
    n_invariant = 0

    with _open(vcf_path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                vcf_samples = header[9:]
                missing = [s for s in sample_list if s not in vcf_samples]
                if missing:
                    raise SystemExit(f"VCF missing samples: {missing}")
                sample_cols = [9 + vcf_samples.index(s) for s in sample_list]
                continue
            if not line or line.startswith("#"):
                continue
            if sample_cols is None:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            alt_field = fields[4]
            if alt_field in (".", "") or "," in alt_field:
                continue
            row = np.empty(n, dtype=np.int8)
            any_missing = False
            for j, idx in enumerate(sample_cols):
                col = fields[idx]
                gt = col if len(col) == 1 else col.split(":", 1)[0]
                if gt == "." or gt == "./." or gt == ".|.":
                    any_missing = True
                    break
                tok = gt.replace("|", "/").split("/")[0]
                if tok == ".":
                    any_missing = True
                    break
                row[j] = 0 if tok == "0" else 1
            if any_missing:
                n_skipped_missing += 1
                continue
            ac = int(row.sum())
            if 0 < ac < n:
                rows.append(row)
            else:
                n_invariant += 1

    G = np.array(rows, dtype=np.int8)
    return G, {"n_skipped_missing": n_skipped_missing, "n_invariant": n_invariant,
               "n_segregating": int(G.shape[0]) if G.size else 0}


def main():
    args = parse_args()
    samples = [s.strip() for s in args.samples.split(",") if s.strip()]
    n = len(samples)
    print(f"Loading {n} samples from VCF...", file=sys.stderr)
    G, stats = build_genotype_matrix(args.vcf, samples)
    print(f"Genotype matrix: {G.shape}  stats={stats}", file=sys.stderr)
    n_sites = G.shape[0]

    # ---- Pairwise IBS ----------------------------------------------------
    ibs = np.eye(n)
    for i, j in combinations(range(n), 2):
        ibs[i, j] = ibs[j, i] = float(np.mean(G[:, i] == G[:, j]))

    ibs_df = pd.DataFrame(ibs, index=samples, columns=samples)
    ibs_df.to_csv(args.ibs_matrix, sep="\t")

    pair_rows = []
    for i, j in combinations(range(n), 2):
        ibs_ij = float(ibs[i, j])
        pair_rows.append({
            "sample_A": samples[i],
            "sample_B": samples[j],
            "IBS": round(ibs_ij, 6),
            "Hamming_distance": round(1.0 - ibs_ij, 6),
            "n_diff_sites": int(np.sum(G[:, i] != G[:, j])),
        })
    pair_df = pd.DataFrame(pair_rows).sort_values("IBS", ascending=False)
    pair_df.to_csv(args.ibs_summary, sep="\t", index=False)

    print("\nPairwise IBS (sorted, top 5):", file=sys.stderr)
    print(pair_df.head(5).to_string(index=False), file=sys.stderr)
    print("\nMean IBS over all pairs:", float(pair_df["IBS"].mean()), file=sys.stderr)
    print("SD IBS:", float(pair_df["IBS"].std()), file=sys.stderr)
    top_pair = pair_df.iloc[0]
    z_top = (top_pair["IBS"] - pair_df["IBS"].mean()) / pair_df["IBS"].std()
    print(f"Top pair: {top_pair['sample_A']} / {top_pair['sample_B']}  "
          f"IBS={top_pair['IBS']:.4f}  z={z_top:.2f}", file=sys.stderr)

    # ---- bin-2 decomposition ---------------------------------------------
    # Folded bin 2 (n=6): allele count == 2 OR allele count == 4. For each
    # such site, count which sample-pair carries the minor allele (AC=2 case)
    # or which pair is REF (AC=4 case, pair is the absent allele = minor).
    ac = G.sum(axis=1)
    bin2_mask = (ac == 2) | (ac == 4)
    n_bin2 = int(bin2_mask.sum())

    pair_counts = {(samples[i], samples[j]): 0 for i, j in combinations(range(n), 2)}
    for site_idx in np.where(bin2_mask)[0]:
        site = G[site_idx]
        ac_here = int(site.sum())
        if ac_here == 2:
            carriers = tuple(sorted(samples[k] for k in np.where(site == 1)[0]))
        else:  # ac_here == 4: the minor allele is REF, carried by the 2 zeros
            carriers = tuple(sorted(samples[k] for k in np.where(site == 0)[0]))
        if carriers in pair_counts:
            pair_counts[carriers] += 1

    bin2_rows = []
    for (a, b), c in pair_counts.items():
        bin2_rows.append({
            "sample_A": a,
            "sample_B": b,
            "bin2_pair_owned_count": c,
            "bin2_pair_owned_fraction": round(c / n_bin2, 4) if n_bin2 else 0.0,
        })
    bin2_df = pd.DataFrame(bin2_rows).sort_values("bin2_pair_owned_count",
                                                  ascending=False)
    bin2_df.to_csv(args.bin2_decomp, sep="\t", index=False)

    print(f"\nbin-2 (folded) total sites: {n_bin2}", file=sys.stderr)
    print("Top pair owners:", file=sys.stderr)
    print(bin2_df.head(5).to_string(index=False), file=sys.stderr)
    expected_uniform = n_bin2 / len(pair_counts)
    print(f"Uniform expectation per pair: {expected_uniform:.1f}", file=sys.stderr)
    top_pair_count = int(bin2_df.iloc[0]["bin2_pair_owned_count"])
    print(f"Top pair: {bin2_df.iloc[0]['sample_A']} / "
          f"{bin2_df.iloc[0]['sample_B']}  count={top_pair_count}  "
          f"~{top_pair_count / expected_uniform:.2f}x uniform",
          file=sys.stderr)

    # ---- Plot ------------------------------------------------------------
    fig, (ax_h, ax_b) = plt.subplots(1, 2, figsize=(14, 5),
                                      gridspec_kw={"width_ratios": [1, 1.2]})

    # Heatmap of IBS (with diagonal as mean(IBS) so colormap focuses on
    # off-diagonal contrast).
    masked = ibs.copy()
    np.fill_diagonal(masked, np.nan)
    im = ax_h.imshow(masked, cmap="viridis_r", aspect="equal")
    ax_h.set_xticks(range(n))
    ax_h.set_yticks(range(n))
    ax_h.set_xticklabels(samples, rotation=45, ha="right")
    ax_h.set_yticklabels(samples)
    ax_h.set_title(f"Pairwise IBS over {n_sites:,} 4D segregating sites\n"
                   f"(diagonal masked; viridis_r → darker = more similar)")
    for i in range(n):
        for j in range(n):
            if i != j:
                ax_h.text(j, i, f"{ibs[i, j]:.3f}", ha="center", va="center",
                          color="white", fontsize=8)
    fig.colorbar(im, ax=ax_h, label="IBS")

    # bar plot of pair-owned bin-2 counts.
    labels = [f"{r['sample_A']}\n{r['sample_B']}" for _, r in bin2_df.iterrows()]
    counts = bin2_df["bin2_pair_owned_count"].values
    ax_b.bar(np.arange(len(labels)), counts, color="steelblue",
             edgecolor="black")
    ax_b.axhline(expected_uniform, color="red", lw=1.5, ls="--",
                 label=f"uniform expectation ({expected_uniform:.1f}/pair)")
    ax_b.set_xticks(np.arange(len(labels)))
    ax_b.set_xticklabels(labels, rotation=80, ha="right", fontsize=7)
    ax_b.set_ylabel("Sites where this pair owns the bin-2 minor allele")
    ax_b.set_title(
        f"bin-2 (folded; AC=2 or AC=4) ownership per pair\n"
        f"total sites in bin-2: {n_bin2}"
    )
    ax_b.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(args.plot, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {args.ibs_matrix}, {args.ibs_summary}, {args.bin2_decomp}, "
          f"{args.plot}", file=sys.stderr)


if __name__ == "__main__":
    main()
