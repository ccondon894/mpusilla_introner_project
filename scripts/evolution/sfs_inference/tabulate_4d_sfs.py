"""Tabulate the folded 1D SFS of 4D sites for Group 1.

Reads the all-sites 4D VCF, restricts to the requested haploid G1 samples,
drops sites with any missing genotype in those samples, and tabulates a
folded SFS of length n+1 (where n = number of samples). Saves the
moments.Spectrum array as .npy and writes a summary TSV-ish text file.
"""

import argparse
import gzip
import sys

import numpy as np
import matplotlib.pyplot as plt
import moments


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf", required=True, help="4D all-sites VCF (bgzipped)")
    p.add_argument(
        "--samples",
        required=True,
        help="Comma-separated list of G1 sample names to include",
    )
    p.add_argument("--output", required=True, help="Path to write folded SFS .npy")
    p.add_argument("--plot", required=True, help="Path to write folded-SFS PNG")
    p.add_argument("--summary", required=True, help="Path to write summary text")
    return p.parse_args()


def _open_vcf(path):
    if path.endswith((".gz", ".bgz")):
        return gzip.open(path, "rt")
    return open(path, "r")


def build_afs(vcf_path, sample_list):
    """Stream-parse the VCF text and build an unfolded 1D allele count spectrum.

    Sequential parse only — no indexing needed. *M. pusilla* is haploid, so
    the sample-column GT is a single allele code (or '.' for missing). We
    accept '|'/'/'/'.' separators defensively and treat any non-zero allele
    as ALT.
    """
    n = len(sample_list)
    afs = np.zeros(n + 1, dtype=np.int64)

    n_total = 0
    n_skipped_multi = 0
    n_skipped_missing = 0
    n_segregating = 0
    n_monomorphic = 0

    sample_cols = None  # indices of requested samples within the data fields

    with _open_vcf(vcf_path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                vcf_samples = header[9:]
                missing = [s for s in sample_list if s not in vcf_samples]
                if missing:
                    raise SystemExit(f"VCF missing requested samples: {missing}")
                # Map each requested sample to its absolute column index.
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
            multi = "," in alt_field

            # GT is the first colon-delimited subfield in each sample column.
            any_missing = False
            alt_count = 0
            for idx in sample_cols:
                col = fields[idx]
                # Fast path: single-character haploid GTs ('0', '1', '.').
                if len(col) == 1:
                    gt = col
                else:
                    gt = col.split(":", 1)[0]
                if gt == "." or gt == "./." or gt == ".|.":
                    any_missing = True
                    break
                # Haploid: gt is a single allele code; diploid would be like
                # '0/0'. Tokenize on '/' and '|' just in case.
                for tok in gt.replace("|", "/").split("/"):
                    if tok == ".":
                        any_missing = True
                        break
                    if tok != "0":
                        alt_count += 1
                if any_missing:
                    break

            if any_missing:
                n_skipped_missing += 1
                continue

            if alt_field in (".", ""):
                # Invariant site — REF only.
                afs[0] += 1
                n_monomorphic += 1
                continue

            if multi:
                n_skipped_multi += 1
                continue

            afs[alt_count] += 1
            if 0 < alt_count < n:
                n_segregating += 1
            else:
                n_monomorphic += 1

    stats = {
        "n_total": n_total,
        "n_skipped_multi": n_skipped_multi,
        "n_skipped_missing": n_skipped_missing,
        "n_segregating": n_segregating,
        "n_monomorphic": n_monomorphic,
        "n_used": int(afs.sum()),
    }
    return afs, stats


def fold_spectrum(afs):
    """Fold via moments and return the folded Spectrum (mask preserved)."""
    fs = moments.Spectrum(afs.astype(float))
    return fs.fold()


def plot_folded(folded, n, plot_path):
    # Moments masks the upper half of folded entries; pull the unmasked indices.
    data = np.ma.filled(folded, 0)
    # Folded SFS only meaningful for i = 1..floor(n/2). Plot those bins only.
    half = n // 2
    bins = np.arange(1, half + 1)
    counts = data[1 : half + 1]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(bins, counts, color="steelblue", edgecolor="black")
    ax.set_xticks(bins)
    ax.set_xlabel("Minor allele count (folded)")
    ax.set_ylabel("Number of 4D sites")
    ax.set_title(f"Folded 4D SFS (Group 1, n={n} haploid)")
    for x, y in zip(bins, counts):
        if y > 0:
            ax.text(x, y, f"{int(y)}", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(plot_path, dpi=200)
    plt.close(fig)


def main():
    args = parse_args()
    sample_list = [s.strip() for s in args.samples.split(",") if s.strip()]
    n = len(sample_list)
    if n < 2:
        raise SystemExit("Need at least 2 samples")

    print(f"Building 4D SFS for {n} samples: {sample_list}", file=sys.stderr)
    afs, stats = build_afs(args.vcf, sample_list)
    folded = fold_spectrum(afs)

    # Save the folded Spectrum's underlying ndarray (with the mask attribute
    # discarded — downstream rebuilds the Spectrum). Also save the unfolded
    # for reproducibility.
    np.savez(
        args.output if args.output.endswith(".npz") else args.output.replace(".npy", ".npz"),
        unfolded=afs,
        folded=np.ma.filled(folded, 0).astype(float),
        folded_mask=folded.mask if hasattr(folded, "mask") else np.zeros_like(afs, dtype=bool),
        n_samples=n,
        sample_names=np.array(sample_list),
    )
    # Also drop a plain .npy of the unfolded spectrum at the requested path
    # (downstream loads it as the canonical input).
    np.save(args.output, afs)

    plot_folded(folded, n, args.plot)

    with open(args.summary, "w") as fh:
        fh.write(f"n_samples\t{n}\n")
        fh.write(f"samples\t{','.join(sample_list)}\n")
        for k, v in stats.items():
            fh.write(f"{k}\t{v}\n")
        fh.write("\n# Unfolded SFS (allele count -> n_sites)\n")
        for i, c in enumerate(afs):
            fh.write(f"unfolded[{i}]\t{int(c)}\n")
        fh.write("\n# Folded SFS (minor allele count -> n_sites)\n")
        data = np.ma.filled(folded, 0)
        for i in range(n // 2 + 1):
            fh.write(f"folded[{i}]\t{int(data[i])}\n")

    print("Stats:", stats, file=sys.stderr)
    print(f"Wrote {args.output}, {args.plot}, {args.summary}", file=sys.stderr)


if __name__ == "__main__":
    main()
