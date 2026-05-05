"""Tabulate the non-introner-intron two-state SFS.

Adapter for the non-introner-intron genotype matrix produced by
`workflow/rules/13_non_introner_introns.smk`. The matrix is in WIDE format
(one row per locus, one column per sample), unlike the introner
`genotype_matrix.final.tsv` which is in long format.

Convention (matches introner pipeline):
    presence == 1  → present (carrier)
    presence == 2  → absent
    presence == 3  → missing

Outputs the same files as `tabulate_introner_sfs.py` so downstream scripts
(fit_gain_loss.py, fit_bin_level_model.py) work without changes.
"""

import argparse
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--matrix", required=True,
                   help="intron_genotype_matrix.tsv (wide format)")
    p.add_argument("--samples", required=True,
                   help="Comma-separated subset of sample names")
    p.add_argument("--output", required=True, help="Output SFS .npy (length n+1)")
    p.add_argument("--summary", required=True, help="Output summary TSV")
    p.add_argument("--plot", required=True, help="Output PNG")
    return p.parse_args()


def main():
    args = parse_args()
    sample_list = [s.strip() for s in args.samples.split(",") if s.strip()]
    n = len(sample_list)

    df = pd.read_csv(args.matrix, sep="\t")
    n_loci_total = len(df)
    print(f"Loaded matrix: {n_loci_total} loci × "
          f"{df.shape[1]} columns", file=sys.stderr)

    missing_cols = [s for s in sample_list if s not in df.columns]
    if missing_cols:
        raise SystemExit(f"Samples missing from matrix: {missing_cols}")

    sub = df[sample_list].copy()

    has_missing = (sub == 3).any(axis=1)
    n_dropped_missing = int(has_missing.sum())
    sub = sub[~has_missing]
    n_loci_called = len(sub)

    present_counts = (sub == 1).sum(axis=1).astype(int)

    sfs_full = np.zeros(n + 1, dtype=np.int64)
    for k in range(n + 1):
        sfs_full[k] = int((present_counts == k).sum())

    fixed_absent = int(sfs_full[0])
    fixed_present = int(sfs_full[n])
    polymorphic = int(sfs_full[1:n].sum())

    np.save(args.output, sfs_full)

    rows = [
        ("n_samples", n),
        ("samples", ",".join(sample_list)),
        ("status_filter", "(none — non-introner intron matrix has no within_group_status column)"),
        ("n_loci_total", n_loci_total),
        ("n_loci_status_filtered", n_loci_total),  # no filter applied
        ("n_loci_dropped_missing", n_dropped_missing),
        ("n_loci_called", n_loci_called),
        ("fixed_present", fixed_present),
        ("fixed_absent", fixed_absent),
        ("polymorphic", polymorphic),
        ("fixed_present_to_absent_ratio",
         f"{fixed_present / fixed_absent:.6f}" if fixed_absent > 0
         else "inf (fixed_absent = 0 by ascertainment — reference catalog only includes "
              "introns present in CCMP1545)"),
    ]
    for i in range(n + 1):
        rows.append((f"sfs[{i}]", int(sfs_full[i])))

    with open(args.summary, "w") as fh:
        fh.write("key\tvalue\n")
        for k, v in rows:
            fh.write(f"{k}\t{v}\n")

    bins = np.arange(1, n)
    counts = sfs_full[1:n]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(bins, counts, color="darkorange", edgecolor="black", alpha=0.75)
    ax.set_xticks(bins)
    ax.set_xlabel("Present-allele count in subset")
    ax.set_ylabel("Number of polymorphic non-introner intron loci")
    ax.set_title(
        f"Non-introner intron SFS (n={n}; "
        f"fixed-P={fixed_present}, fixed-A={fixed_absent})"
    )
    for x, y in zip(bins, counts):
        if y > 0:
            ax.text(x, y, str(int(y)), ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(args.plot, dpi=200)
    plt.close(fig)

    print(
        f"n_loci_total={n_loci_total}  called={n_loci_called}  "
        f"dropped_missing={n_dropped_missing}",
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
