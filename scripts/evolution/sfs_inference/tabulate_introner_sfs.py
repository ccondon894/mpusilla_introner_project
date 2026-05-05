"""Tabulate the introner two-state SFS for Group 1.

Reads `genotype_matrix.final.tsv`, restricts to G1 samples, drops loci with
any missing G1 call, applies the within_group_status filter, and tabulates
counts of present-allele frequency 0..n.

Convention (per docs/haplotype_diversity_diagnostics.md):
    presence == 1  → present (carrier)
    presence == 2  → absent
    presence == 3  → missing
"""

import argparse
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--matrix", required=True, help="genotype_matrix.final.tsv")
    p.add_argument("--samples", required=True, help="Comma-separated G1 sample list")
    p.add_argument("--output", required=True, help="Output SFS .npy (length n+1)")
    p.add_argument("--summary", required=True, help="Output summary TSV")
    p.add_argument("--plot", required=True, help="Output PNG")
    p.add_argument(
        "--status_filter",
        default="consistent,singleton",
        help="Comma-separated within_group_status values to keep",
    )
    return p.parse_args()


def main():
    args = parse_args()
    sample_list = [s.strip() for s in args.samples.split(",") if s.strip()]
    n = len(sample_list)
    keep_status = set(s.strip() for s in args.status_filter.split(",") if s.strip())

    df = pd.read_csv(args.matrix, sep="\t")
    print(f"Loaded matrix: {len(df)} rows", file=sys.stderr)

    df = df[df["sample"].isin(sample_list)].copy()
    print(f"After G1 sample filter: {len(df)} rows", file=sys.stderr)

    # Pivot: one row per ortholog, one column per sample, presence as values.
    pivot = df.pivot(index="ortholog_id", columns="sample", values="presence")
    # Reorder columns to match the requested sample list. If any are missing
    # entirely, pandas leaves them out — surface that.
    missing_cols = [s for s in sample_list if s not in pivot.columns]
    if missing_cols:
        raise SystemExit(f"Samples missing from matrix: {missing_cols}")
    pivot = pivot[sample_list]

    # Per-locus within_group_status: each ortholog should have a single value
    # across its sample rows. Take the first non-null.
    status = (
        df.groupby("ortholog_id")["within_group_status"].first().reindex(pivot.index)
    )

    n_loci_total = len(pivot)

    # Apply within_group_status filter.
    keep = status.isin(keep_status)
    pivot_kept = pivot[keep]
    status_kept = status[keep]
    n_loci_status = len(pivot_kept)

    # Drop loci with any missing G1 (presence == 3).
    has_missing = (pivot_kept == 3).any(axis=1)
    pivot_called = pivot_kept[~has_missing]
    n_loci_called = len(pivot_called)
    n_dropped_missing = int(has_missing.sum())

    # Per-locus present count (presence == 1).
    present_counts = (pivot_called == 1).sum(axis=1).astype(int)

    sfs_full = np.zeros(n + 1, dtype=np.int64)
    for k in range(n + 1):
        sfs_full[k] = int((present_counts == k).sum())

    fixed_absent = int(sfs_full[0])
    fixed_present = int(sfs_full[n])
    polymorphic = int(sfs_full[1:n].sum())

    # Save the full-length SFS (includes fixed bins). Downstream uses
    # sfs_full[1:n] for the segregating-only fit.
    np.save(args.output, sfs_full)

    # Summary TSV (key/value).
    rows = [
        ("n_samples", n),
        ("samples", ",".join(sample_list)),
        ("status_filter", ",".join(sorted(keep_status))),
        ("n_loci_total", n_loci_total),
        ("n_loci_status_filtered", n_loci_status),
        ("n_loci_dropped_missing", n_dropped_missing),
        ("n_loci_called", n_loci_called),
        ("fixed_present", fixed_present),
        ("fixed_absent", fixed_absent),
        ("polymorphic", polymorphic),
        ("fixed_present_to_absent_ratio",
         f"{fixed_present / fixed_absent:.6f}" if fixed_absent > 0 else "inf"),
    ]
    for i in range(n + 1):
        rows.append((f"sfs[{i}]", int(sfs_full[i])))

    with open(args.summary, "w") as fh:
        fh.write("key\tvalue\n")
        for k, v in rows:
            fh.write(f"{k}\t{v}\n")

    # Plot the segregating SFS (i = 1..n-1).
    bins = np.arange(1, n)
    counts = sfs_full[1:n]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(bins, counts, color="purple", edgecolor="black", alpha=0.75)
    ax.set_xticks(bins)
    ax.set_xlabel("Present-allele count in G1")
    ax.set_ylabel("Number of polymorphic introner loci")
    ax.set_title(
        f"Introner two-state SFS (n={n}; fixed-P={fixed_present}, fixed-A={fixed_absent})"
    )
    for x, y in zip(bins, counts):
        if y > 0:
            ax.text(x, y, str(int(y)), ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(args.plot, dpi=200)
    plt.close(fig)

    print(
        f"n_loci_total={n_loci_total}  status_kept={n_loci_status}  "
        f"called={n_loci_called}  dropped_missing={n_dropped_missing}",
        file=sys.stderr,
    )
    print(
        f"fixed_P={fixed_present}  fixed_A={fixed_absent}  poly={polymorphic}",
        file=sys.stderr,
    )
    print(f"SFS (i=0..{n}): {sfs_full.tolist()}", file=sys.stderr)
    print(f"Wrote {args.output}, {args.summary}, {args.plot}", file=sys.stderr)


if __name__ == "__main__":
    main()
