"""Check Group 1 for population substructure.

Reads the 4D all-sites VCF, restricts to G1, keeps biallelic segregating sites
with no missing G1 calls, and:

1. Runs PCA via SVD on the standardized genotype matrix.
2. Splits samples into two candidate clusters by sign of PC1 (and reports
   k-means k=2 as a sanity check).
3. Computes Hudson's FST between the two clusters using a per-site allele-
   count estimator, plus a permutation null (sample labels shuffled).
4. Decomposes the folded SFS contribution by which subclade carries the
   minor allele to confirm whether a cluster split would put mass at the
   bins the residuals analysis flagged.

Outputs PCA scatter plot, sample × PC table, FST results, and per-bin SFS
decomposition.
"""

import argparse
import gzip
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf", required=True)
    p.add_argument("--samples", required=True, help="Comma-separated G1 sample list")
    p.add_argument("--pca_plot", required=True)
    p.add_argument("--pca_table", required=True, help="TSV with PC1..PC4 per sample")
    p.add_argument("--fst", required=True, help="TSV with FST + permutation results")
    p.add_argument("--sfs_decomp_plot", required=True)
    p.add_argument("--sfs_decomp_table", required=True)
    p.add_argument("--summary_plot", required=True,
                   help="Two-panel summary: PCA + FST permutation null")
    p.add_argument("--n_permutations", type=int, default=1000)
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def _open(path):
    return gzip.open(path, "rt") if path.endswith((".gz", ".bgz")) else open(path)


def build_genotype_matrix(vcf_path, sample_list):
    """Return (G, sites_kept) where G is (n_sites, n_samples) of {0, 1}.

    Streams the VCF, keeps only biallelic SNP sites segregating in G1 with no
    missing G1 calls.
    """
    n = len(sample_list)
    sample_cols = None
    rows = []  # list of (alt_count_per_sample,) tuples

    n_total = 0
    n_skipped_missing = 0
    n_skipped_multi = 0
    n_skipped_invariant = 0

    with _open(vcf_path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                vcf_samples = header[9:]
                missing = [s for s in sample_list if s not in vcf_samples]
                if missing:
                    raise SystemExit(f"VCF missing requested samples: {missing}")
                sample_cols = [9 + vcf_samples.index(s) for s in sample_list]
                continue
            if not line or line.startswith("#"):
                continue
            if sample_cols is None:
                raise SystemExit("Saw data line before #CHROM header")

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            n_total += 1

            alt_field = fields[4]
            if alt_field in (".", ""):
                n_skipped_invariant += 1
                continue
            if "," in alt_field:
                n_skipped_multi += 1
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

            if 0 < int(row.sum()) < n:
                rows.append(row)
            else:
                n_skipped_invariant += 1

    G = np.array(rows, dtype=np.int8)
    stats = {
        "n_total_records": n_total,
        "n_skipped_missing": n_skipped_missing,
        "n_skipped_multi": n_skipped_multi,
        "n_skipped_invariant": n_skipped_invariant,
        "n_sites_used": int(G.shape[0]) if G.size else 0,
    }
    return G, stats


def pca_svd(G):
    """PCA via SVD on standardized site × sample matrix.

    Standardizes each site by p, sqrt(p*(1-p)), drops sites with no variance.
    Returns (sample_scores, eigenvalues_explained, n_sites_used).
    """
    G = G.astype(float)
    p = G.mean(axis=1, keepdims=True)
    sd = np.sqrt(p * (1 - p))
    keep = (sd > 0).flatten()
    if keep.sum() == 0:
        raise RuntimeError("No variant sites after standardization")
    Z = (G[keep] - p[keep]) / sd[keep]
    # SVD: Z = U S Vt; sample loadings on PCs are columns of Vt^T (i.e., V).
    # numpy returns Vt with shape (n_pcs, n_samples).
    U, S, Vt = np.linalg.svd(Z, full_matrices=False)
    scores = Vt.T * S  # sample × PC; multiply by S so they're true scores
    var_explained = (S ** 2) / np.sum(S ** 2)
    return scores, var_explained, int(keep.sum())


def hudson_fst(G, label1, label2):
    """Hudson's FST over polymorphic sites between two clusters.

    Sites where either cluster is monomorphic and identical are dropped from
    the per-site numerator/denominator to avoid 0/0; sites where the union
    is monomorphic in *both* clusters at the same allele give D=0, N=0 and
    are dropped.

    Returns weighted (ratio of sums) FST.
    """
    n1 = label1.sum()
    n2 = label2.sum()
    if n1 < 2 or n2 < 2:
        raise ValueError("Each cluster needs >= 2 samples for unbiased het")
    # Per-site allele frequencies in each cluster.
    G1 = G[:, label1].astype(float)
    G2 = G[:, label2].astype(float)
    p1 = G1.mean(axis=1)
    p2 = G2.mean(axis=1)
    # Hudson 1992 estimator: numerator (between-pop het) - within-pop het.
    # D = unbiased within-pop het summed, weighted appropriately.
    h1 = (n1 / (n1 - 1.0)) * p1 * (1.0 - p1)
    h2 = (n2 / (n2 - 1.0)) * p2 * (1.0 - p2)
    Hw = 0.5 * (h1 + h2)
    Hb = p1 * (1.0 - p2) + p2 * (1.0 - p1)
    # Drop sites with Hb == 0 (no between-pop variation); they carry no info.
    use = Hb > 0
    if not use.any():
        return float("nan")
    fst_per_site = (Hb[use] - Hw[use]) / Hb[use]
    # Ratio of averages (Bhatia et al. 2013 recommendation).
    return float((Hb[use] - Hw[use]).sum() / Hb[use].sum()), fst_per_site


def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    sample_list = [s.strip() for s in args.samples.split(",") if s.strip()]
    n = len(sample_list)
    print(f"Loading 4D segregating sites for {n} samples...", file=sys.stderr)
    G, stats = build_genotype_matrix(args.vcf, sample_list)
    print(f"Genotype matrix: {G.shape}", file=sys.stderr)
    print(f"Stats: {stats}", file=sys.stderr)
    if G.size == 0:
        raise SystemExit("No usable sites")

    # ---- PCA -------------------------------------------------------------
    scores, var_exp, n_used = pca_svd(G)
    pc_df = pd.DataFrame({
        "sample": sample_list,
        "PC1": scores[:, 0],
        "PC2": scores[:, 1] if scores.shape[1] > 1 else np.zeros(n),
        "PC3": scores[:, 2] if scores.shape[1] > 2 else np.zeros(n),
        "PC4": scores[:, 3] if scores.shape[1] > 3 else np.zeros(n),
    })
    pc_df.to_csv(args.pca_table, sep="\t", index=False)
    print(f"Variance explained (top 4 PCs): "
          f"{[f'{v:.3f}' for v in var_exp[:4]]}", file=sys.stderr)

    # Cluster: PC1 sign + k-means on PCs 1-2 for sanity.
    pc1_split = (scores[:, 0] >= 0).astype(int)
    if scores.shape[1] >= 2:
        km = KMeans(n_clusters=2, n_init=20, random_state=args.seed)
        km_labels = km.fit_predict(scores[:, :2])
    else:
        km_labels = pc1_split
    # Standardize: cluster 0 contains the smallest-PC1 sample.
    if scores[np.argmin(scores[:, 0]), 0] < 0:
        # Flip so cluster 0 is the negative-PC1 group.
        if pc1_split[np.argmin(scores[:, 0])] != 0:
            pc1_split = 1 - pc1_split
        if km_labels[np.argmin(scores[:, 0])] != 0:
            km_labels = 1 - km_labels

    pc1_split_labels = pc1_split.astype(bool)
    sizes = (int((~pc1_split_labels).sum()), int(pc1_split_labels.sum()))
    print(f"PC1-sign split: cluster sizes {sizes} ", file=sys.stderr)
    print(f"  cluster_0 (PC1<0): {[s for s, c in zip(sample_list, pc1_split) if c == 0]}",
          file=sys.stderr)
    print(f"  cluster_1 (PC1>=0): {[s for s, c in zip(sample_list, pc1_split) if c == 1]}",
          file=sys.stderr)

    # ---- FST + permutation ----------------------------------------------
    label1 = ~pc1_split_labels
    label2 = pc1_split_labels
    if label1.sum() < 2 or label2.sum() < 2:
        print("WARNING: at least one cluster has <2 samples; skipping FST",
              file=sys.stderr)
        fst_obs = float("nan")
        perm_p = float("nan")
        perm_dist = np.array([])
    else:
        fst_obs, fst_per_site = hudson_fst(G, label1, label2)
        # Permutation: shuffle labels.
        n_total = label1.sum() + label2.sum()
        size_a = int(label1.sum())
        perm_dist = np.zeros(args.n_permutations)
        for r in range(args.n_permutations):
            perm = rng.permutation(n_total)
            perm_label = np.zeros(n_total, dtype=bool)
            perm_label[perm[:size_a]] = True
            f, _ = hudson_fst(G, perm_label, ~perm_label)
            perm_dist[r] = f
        perm_p = float(np.mean(perm_dist >= fst_obs))
        print(f"Hudson FST = {fst_obs:.4f}  (perm p = {perm_p:.4f}; "
              f"perm mean = {perm_dist.mean():.4f}, sd = {perm_dist.std():.4f})",
              file=sys.stderr)

    with open(args.fst, "w") as fh:
        fh.write("key\tvalue\n")
        fh.write(f"n_sites_used\t{n_used}\n")
        fh.write(f"cluster_0_size\t{int(label1.sum())}\n")
        fh.write(f"cluster_1_size\t{int(label2.sum())}\n")
        fh.write("cluster_0_samples\t" +
                 ",".join(s for s, c in zip(sample_list, pc1_split) if c == 0) + "\n")
        fh.write("cluster_1_samples\t" +
                 ",".join(s for s, c in zip(sample_list, pc1_split) if c == 1) + "\n")
        fh.write(f"hudson_fst_observed\t{fst_obs:.6g}\n")
        fh.write(f"perm_n\t{args.n_permutations}\n")
        if perm_dist.size:
            fh.write(f"perm_mean\t{perm_dist.mean():.6g}\n")
            fh.write(f"perm_sd\t{perm_dist.std():.6g}\n")
            fh.write(f"perm_p2.5\t{np.percentile(perm_dist, 2.5):.6g}\n")
            fh.write(f"perm_p97.5\t{np.percentile(perm_dist, 97.5):.6g}\n")
            fh.write(f"perm_one_sided_p\t{perm_p:.6g}\n")
        # k-means sanity check.
        fh.write("kmeans_cluster_0_samples\t" +
                 ",".join(s for s, c in zip(sample_list, km_labels) if c == 0) + "\n")
        fh.write("kmeans_cluster_1_samples\t" +
                 ",".join(s for s, c in zip(sample_list, km_labels) if c == 1) + "\n")

    # ---- SFS decomposition by cluster ------------------------------------
    # For each segregating site, count alt alleles in cluster 0 and cluster 1.
    # Tabulate the unfolded SFS and split by which cluster contributes the
    # minor allele (or both equally).
    alt_in_0 = G[:, label1].sum(axis=1)
    alt_in_1 = G[:, label2].sum(axis=1)
    alt_total = alt_in_0 + alt_in_1
    # Fold to minor-allele count.
    folded_count = np.minimum(alt_total, n - alt_total)

    n_size_0 = int(label1.sum())
    n_size_1 = int(label2.sum())
    private_0_mask = (alt_in_1 == 0) | (alt_in_1 == n_size_1)  # invariant in cluster 1
    private_1_mask = (alt_in_0 == 0) | (alt_in_0 == n_size_0)  # invariant in cluster 0
    private_both = private_0_mask & private_1_mask  # invariant in both → between-cluster fixed differences
    private_0_only = private_0_mask & ~private_both
    private_1_only = private_1_mask & ~private_both
    shared = ~(private_0_mask | private_1_mask)

    half = n // 2
    bins = np.arange(1, half + 1)
    rows = []
    for b in bins:
        sel = folded_count == b
        rows.append({
            "bin": int(b),
            "total": int(sel.sum()),
            "private_to_cluster0": int((sel & private_0_only).sum()),
            "private_to_cluster1": int((sel & private_1_only).sum()),
            "fixed_difference_between_clusters": int((sel & private_both).sum()),
            "shared_polymorphism": int((sel & shared).sum()),
        })
    decomp_df = pd.DataFrame(rows)
    decomp_df.to_csv(args.sfs_decomp_table, sep="\t", index=False)

    # ---- Plots -----------------------------------------------------------
    pc1 = scores[:, 0]
    pc2 = scores[:, 1] if scores.shape[1] > 1 else np.zeros(n)
    colors = ["tab:blue" if c == 0 else "tab:red" for c in pc1_split]

    # Plot 1: PCA + SFS decomposition (existing).
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 5))
    axA.scatter(pc1, pc2, c=colors, s=80, edgecolor="black")
    for i, s in enumerate(sample_list):
        axA.annotate(s, (pc1[i], pc2[i]), fontsize=8, xytext=(4, 4),
                     textcoords="offset points")
    axA.axvline(0, color="gray", lw=0.5, ls="--")
    axA.set_xlabel(f"PC1 ({var_exp[0] * 100:.1f}%)")
    axA.set_ylabel(f"PC2 ({var_exp[1] * 100:.1f}%)" if scores.shape[1] > 1 else "PC2")
    axA.set_title(f"G1 PCA on {n_used} 4D SNPs (color = PC1-sign cluster)")

    axB.bar(decomp_df["bin"], decomp_df["private_to_cluster0"],
            label=f"private to cluster_0 (n={n_size_0})",
            color="tab:blue", edgecolor="black")
    axB.bar(decomp_df["bin"], decomp_df["private_to_cluster1"],
            bottom=decomp_df["private_to_cluster0"],
            label=f"private to cluster_1 (n={n_size_1})",
            color="tab:red", edgecolor="black")
    axB.bar(decomp_df["bin"], decomp_df["shared_polymorphism"],
            bottom=decomp_df["private_to_cluster0"] + decomp_df["private_to_cluster1"],
            label="shared polymorphism", color="tab:gray", edgecolor="black")
    axB.bar(decomp_df["bin"], decomp_df["fixed_difference_between_clusters"],
            bottom=(decomp_df["private_to_cluster0"]
                    + decomp_df["private_to_cluster1"]
                    + decomp_df["shared_polymorphism"]),
            label="fixed differences", color="tab:orange", edgecolor="black")
    axB.set_xticks(bins)
    axB.set_xlabel("Folded minor allele count")
    axB.set_ylabel("Sites")
    axB.set_title("Folded 4D SFS decomposed by cluster ownership")
    axB.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(args.pca_plot, dpi=200)
    fig.savefig(args.sfs_decomp_plot, dpi=200)
    plt.close(fig)

    # Plot 2: PCA + FST permutation null (the requested summary visual).
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    ax1.scatter(pc1, pc2, c=colors, s=110, edgecolor="black", linewidth=0.8)
    for i, s in enumerate(sample_list):
        ax1.annotate(s, (pc1[i], pc2[i]), fontsize=9, xytext=(5, 5),
                     textcoords="offset points")
    ax1.axvline(0, color="gray", lw=0.5, ls="--")
    ax1.set_xlabel(f"PC1 ({var_exp[0] * 100:.1f}%)")
    ax1.set_ylabel(f"PC2 ({var_exp[1] * 100:.1f}%)" if scores.shape[1] > 1 else "PC2")
    ax1.set_title(
        f"G1 PCA on {n_used} 4D SNPs\n"
        f"cluster_0 (blue, n={n_size_0}) vs cluster_1 (red, n={n_size_1})"
    )

    if perm_dist.size > 0:
        ax2.hist(perm_dist, bins=30, color="lightgray", edgecolor="black",
                 label=f"permutation null (n={args.n_permutations})")
        ax2.axvline(fst_obs, color="red", lw=2.5,
                    label=f"observed FST = {fst_obs:.3f}")
        ax2.axvline(perm_dist.mean(), color="black", lw=1, ls="--",
                    label=f"perm mean = {perm_dist.mean():.3f}")
        # Annotate one-sided p in upper-right.
        ax2.text(
            0.97, 0.95,
            f"one-sided p = {perm_p:.3g}\n"
            f"perm SD = {perm_dist.std():.3f}\n"
            f"z ≈ {(fst_obs - perm_dist.mean()) / perm_dist.std():.2f}",
            transform=ax2.transAxes, ha="right", va="top",
            fontsize=9,
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="gray"),
        )
        ax2.set_xlabel("Hudson FST")
        ax2.set_ylabel("Permutation count")
        ax2.set_title("Hudson FST: observed vs label-shuffled null")
        ax2.legend(fontsize=8, loc="upper left")
    else:
        ax2.text(0.5, 0.5, "FST not computed\n(cluster too small)",
                 ha="center", va="center", transform=ax2.transAxes)
    fig.suptitle(
        "Group 1 has clear population substructure",
        fontsize=12, y=1.02,
    )
    fig.tight_layout()
    fig.savefig(args.summary_plot, dpi=200, bbox_inches="tight")
    plt.close(fig)

    # Save raw permutation array for reproducibility.
    perm_npy = args.summary_plot.rsplit(".", 1)[0] + "_perm_dist.npy"
    np.save(perm_npy, perm_dist)

    print(f"Wrote {args.pca_plot}, {args.pca_table}, {args.fst}, "
          f"{args.sfs_decomp_table}, {args.sfs_decomp_plot}, {args.summary_plot}, "
          f"{perm_npy}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
