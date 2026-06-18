#!/usr/bin/env python3
"""Build RCC1749 pre-insertion-site sequence features for random forest tests."""

import argparse
import bisect
import random
import re
from collections import Counter, defaultdict
from itertools import product
from pathlib import Path

import pandas as pd
from Bio import SeqIO


DEFAULT_WINDOWS = [25, 50, 100]
KMERS_3 = ["".join(p) for p in product("ACGT", repeat=3)]
VALID_BASES = set("ACGT")
REVCOMP_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--genotype-matrix",
        default="../../results/genotyping/genotype_matrix.final.tsv",
    )
    parser.add_argument(
        "--gtf",
        default="../../results/annotations/RCC1749.gtf",
    )
    parser.add_argument(
        "--fa",
        default="../../results/assemblies/RCC1749.vg_paths.fa",
    )
    parser.add_argument("--sample", default="RCC1749")
    parser.add_argument("--reference-sample", default="CCMP1545")
    parser.add_argument("--windows", default="25,50,100")
    parser.add_argument("--controls-per-positive", type=int, default=1)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument(
        "--junction-mask-bp",
        type=int,
        default=3,
        help="Mask this many bases on either side of each candidate insertion midpoint.",
    )
    parser.add_argument(
        "--exclude-contig",
        action="append",
        default=[],
        help="Contig to exclude from positives and controls; can be provided more than once.",
    )
    parser.add_argument(
        "--mummer-coords",
        default=None,
        help=(
            "Optional MUMmer show-coords file. Query contigs with mostly reverse "
            "alignments will have left/right flanks reverse-complemented and swapped."
        ),
    )
    parser.add_argument(
        "--orientation-threshold",
        type=float,
        default=0.8,
        help="Minimum reverse-aligned fraction for treating a query contig as reversed.",
    )
    parser.add_argument(
        "--keep-positives-outside-exons",
        action="store_true",
        help="Keep positive absent sites even if the max window is not fully inside one exon.",
    )
    parser.add_argument(
        "--max-absent-span",
        type=int,
        default=None,
        help="Optional maximum RCC1749 absent interval span to retain as a positive.",
    )
    parser.add_argument(
        "--output",
        default="rf_results/rcc1749_preinsertion_sequence_features.tsv",
    )
    parser.add_argument(
        "--summary",
        default="rf_results/rcc1749_preinsertion_sequence_features.summary.tsv",
    )
    return parser.parse_args()


def parse_windows(value):
    return [int(part.strip()) for part in value.split(",") if part.strip()]


def parse_gtf_attrs(attr_text):
    return dict(re.findall(r'(\S+) "([^"]+)"', attr_text))


def read_exons(gtf_path):
    exons = defaultdict(list)
    with open(gtf_path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            attrs = parse_gtf_attrs(fields[8])
            exons[fields[0]].append({
                "contig": fields[0],
                "start": int(fields[3]) - 1,
                "end": int(fields[4]),
                "gene": attrs.get("gene_id", ""),
                "transcript_id": attrs.get("transcript_id", ""),
                "strand": fields[6],
            })
    for contig in exons:
        exons[contig].sort(key=lambda row: (row["start"], row["end"], row["gene"]))
    return exons


def read_fai(fasta_path):
    sizes = {}
    with open(f"{fasta_path}.fai") as handle:
        for line in handle:
            contig, size, *_ = line.rstrip("\n").split("\t")
            sizes[contig] = int(size)
    return sizes


def build_interval_index(intervals_by_contig):
    index = {}
    for contig, intervals in intervals_by_contig.items():
        merged = []
        for start, end in sorted(intervals):
            if not merged or start > merged[-1][1]:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)
        starts = [start for start, _ in merged]
        index[contig] = (starts, merged)
    return index


def overlaps_interval(index, contig, start, end):
    starts, intervals = index.get(contig, ([], []))
    pos = bisect.bisect_right(starts, start)
    for idx in (pos - 1, pos):
        if 0 <= idx < len(intervals):
            interval_start, interval_end = intervals[idx]
            if start < interval_end and end > interval_start:
                return True
    return False


def containing_exon(exons_by_contig, contig, start, end):
    for exon in exons_by_contig.get(contig, []):
        if exon["start"] <= start and end <= exon["end"]:
            return exon
        if exon["start"] > start:
            break
    return None


def fetch_interval(fasta_index, contig, start, end):
    if contig not in fasta_index:
        return ""
    start = max(0, int(start))
    end = min(int(end), len(fasta_index[contig]))
    if end <= start:
        return ""
    return str(fasta_index[contig].seq[start:end].upper())


def mask_junction_bases(left_flank: str, right_flank: str, n: int = 3) -> tuple[str, str]:
    if n <= 0:
        return left_flank, right_flank
    if len(left_flank) >= n:
        left_flank = left_flank[:-n] + "N" * n
    if len(right_flank) >= n:
        right_flank = "N" * n + right_flank[n:]
    return left_flank, right_flank


def reverse_complement(seq):
    return seq.translate(REVCOMP_TABLE)[::-1].upper()


def read_reverse_oriented_contigs(coords_path, threshold):
    if not coords_path:
        return set()

    by_contig = {}
    with open(coords_path) as handle:
        for line in handle:
            line = line.strip()
            if (
                not line
                or line.startswith("/")
                or line.startswith("NUCMER")
                or line.startswith("[")
            ):
                continue
            parts = line.split("\t")
            if len(parts) < 13:
                continue
            try:
                s2, e2 = int(parts[2]), int(parts[3])
                len1, len2 = int(parts[4]), int(parts[5])
            except ValueError:
                continue
            query_contig = parts[12]
            aligned_bp = min(len1, len2)
            if aligned_bp <= 0:
                continue
            orientation = by_contig.setdefault(query_contig, {"forward": 0, "reverse": 0})
            if s2 > e2:
                orientation["reverse"] += aligned_bp
            else:
                orientation["forward"] += aligned_bp

    reverse_contigs = set()
    for contig, counts in by_contig.items():
        total = counts["forward"] + counts["reverse"]
        if total and counts["reverse"] / total >= threshold:
            reverse_contigs.add(contig)
    return reverse_contigs


def count_overlapping(seq, motif):
    k = len(motif)
    return sum(seq[i:i + k] == motif for i in range(len(seq) - k + 1))


def kmer_frequencies(seq, k=3, kmers=KMERS_3):
    seq = seq.upper()
    counts = Counter()
    total = 0
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if set(kmer) <= VALID_BASES:
            counts[kmer] += 1
            total += 1
    if total == 0:
        return {kmer: 0.0 for kmer in kmers}
    return {kmer: counts[kmer] / total for kmer in kmers}


def longest_homopolymer_run(seq):
    longest = 0
    current = 0
    previous = None
    for base in seq.upper():
        if base == previous and base in VALID_BASES:
            current += 1
        elif base in VALID_BASES:
            current = 1
            previous = base
        else:
            current = 0
            previous = None
        longest = max(longest, current)
    return longest


def scalar_features(seq):
    seq = seq.upper()
    n = len(seq)
    if n == 0:
        return {
            "gc": float("nan"),
            "cpg_density": float("nan"),
            "gt_density": float("nan"),
            "ag_density": float("nan"),
            "longest_homopolymer_run": float("nan"),
        }
    denom_2mer = max(n - 1, 1)
    return {
        "gc": (seq.count("G") + seq.count("C")) / n,
        "cpg_density": count_overlapping(seq, "CG") / denom_2mer,
        "gt_density": count_overlapping(seq, "GT") / denom_2mer,
        "ag_density": count_overlapping(seq, "AG") / denom_2mer,
        "longest_homopolymer_run": longest_homopolymer_run(seq),
    }


def add_window_features(row, fasta_index, windows, junction_mask_bp=0):
    contig = row["contig"]
    site = int(row["site"])
    reverse_oriented = bool(row.get("reverse_oriented_contig", 0))
    features = {}
    for window in windows:
        if reverse_oriented:
            left = reverse_complement(fetch_interval(fasta_index, contig, site, site + window))
            right = reverse_complement(fetch_interval(fasta_index, contig, site - window, site))
        else:
            left = fetch_interval(fasta_index, contig, site - window, site)
            right = fetch_interval(fasta_index, contig, site, site + window)
        left, right = mask_junction_bases(left, right, n=junction_mask_bp)
        combined = left + right

        left_scalars = scalar_features(left)
        right_scalars = scalar_features(right)
        combined_scalars = scalar_features(combined)
        for name, value in left_scalars.items():
            features[f"{name}_left_{window}"] = value
        for name, value in right_scalars.items():
            features[f"{name}_right_{window}"] = value
        for name, value in combined_scalars.items():
            features[f"{name}_combined_{window}"] = value
        for name in left_scalars:
            features[f"{name}_delta_{window}"] = left_scalars[name] - right_scalars[name]

        left_kmers = kmer_frequencies(left)
        right_kmers = kmer_frequencies(right)
        for kmer in KMERS_3:
            features[f"kmer_{kmer}_left_{window}"] = left_kmers[kmer]
            features[f"kmer_{kmer}_right_{window}"] = right_kmers[kmer]
            features[f"kmer_{kmer}_delta_{window}"] = left_kmers[kmer] - right_kmers[kmer]
    return features


def positive_rows(genotype_df, exons_by_contig, fasta_sizes, max_window, args):
    subset = genotype_df.loc[
        genotype_df["sample"].isin([args.reference_sample, args.sample]),
        ["ortholog_id", "sample", "presence"],
    ]
    wide = subset.pivot_table(
        index="ortholog_id",
        columns="sample",
        values="presence",
        aggfunc="first",
    )
    candidate_ids = set(
        wide.loc[
            wide[args.reference_sample].eq(1) & wide[args.sample].eq(2)
        ].index
    )

    rows = genotype_df.loc[
        genotype_df["sample"].eq(args.sample)
        & genotype_df["ortholog_id"].isin(candidate_ids)
    ].copy()
    if args.exclude_contig:
        rows = rows.loc[~rows["contig"].isin(set(args.exclude_contig))].copy()
    rows["span"] = rows["end"].astype(int) - rows["start"].astype(int)
    if args.max_absent_span is not None:
        rows = rows.loc[rows["span"] <= args.max_absent_span].copy()
    rows["site"] = ((rows["start"].astype(int) + rows["end"].astype(int)) // 2).astype(int)

    out_rows = []
    for row in rows.to_dict("records"):
        contig = row["contig"]
        site = int(row["site"])
        if contig not in fasta_sizes or site - max_window < 0 or site + max_window > fasta_sizes[contig]:
            continue
        exon = containing_exon(exons_by_contig, contig, site - max_window, site + max_window)
        if exon is None and not args.keep_positives_outside_exons:
            continue
        out_rows.append({
            "sequence_id": f"{contig}:{site}-{site}:preinsertion:{row['ortholog_id']}",
            "ortholog_id": row["ortholog_id"],
            "sample": args.sample,
            "source": "rcc1749_absent_introner_site",
            "contig": contig,
            "start": int(row["start"]),
            "end": int(row["end"]),
            "site": site,
            "feature_start": site,
            "feature_end": site,
            "gene": exon["gene"] if exon else "",
            "transcript_id": exon["transcript_id"] if exon else "",
            "strand": exon["strand"] if exon else "",
            "family": row.get("family", ""),
            "absent_span": int(row["span"]),
            "exon_start": exon["start"] if exon else pd.NA,
            "exon_end": exon["end"] if exon else pd.NA,
            "reverse_oriented_contig": int(contig in args.reverse_contigs),
            "label": 1,
        })
    return pd.DataFrame(out_rows)


def build_exclusion_index(genotype_df, sample):
    intervals = defaultdict(list)
    rows = genotype_df.loc[genotype_df["sample"].eq(sample)]
    for row in rows.itertuples(index=False):
        intervals[row.contig].append((int(row.start), int(row.end)))
    return build_interval_index(intervals)


def sample_control_rows(
    exons_by_contig,
    fasta_sizes,
    exclusion_index,
    n_controls,
    max_window,
    sample,
    random_state,
    reverse_contigs=None,
    excluded_contigs=None,
):
    rng = random.Random(random_state)
    reverse_contigs = reverse_contigs or set()
    excluded_contigs = excluded_contigs or set()
    candidate_exons = []
    weights = []
    for contig, exons in exons_by_contig.items():
        if contig in excluded_contigs:
            continue
        if contig not in fasta_sizes:
            continue
        for exon in exons:
            low = max(exon["start"] + max_window, 0)
            high = min(exon["end"] - max_window, fasta_sizes[contig])
            if high <= low:
                continue
            candidate_exons.append((exon, low, high))
            weights.append(high - low)

    controls = []
    used_sites = set()
    max_attempts = max(n_controls * 200, 10000)
    for attempt in range(max_attempts):
        if len(controls) >= n_controls:
            break
        exon, low, high = rng.choices(candidate_exons, weights=weights, k=1)[0]
        site = rng.randrange(low, high)
        site_key = (exon["contig"], site)
        if site_key in used_sites:
            continue
        if overlaps_interval(exclusion_index, exon["contig"], site - max_window, site + max_window):
            continue
        used_sites.add(site_key)
        control_id = f"control_{len(controls) + 1:06d}"
        controls.append({
            "sequence_id": f"{exon['contig']}:{site}-{site}:exonic_control:{control_id}",
            "ortholog_id": "",
            "sample": sample,
            "source": "rcc1749_exonic_control",
            "contig": exon["contig"],
            "start": site,
            "end": site,
            "site": site,
            "feature_start": site,
            "feature_end": site,
            "gene": exon["gene"],
            "transcript_id": exon["transcript_id"],
            "strand": exon["strand"],
            "family": "exonic_control",
            "absent_span": pd.NA,
            "exon_start": exon["start"],
            "exon_end": exon["end"],
            "reverse_oriented_contig": int(exon["contig"] in reverse_contigs),
            "label": 0,
        })
    if len(controls) < n_controls:
        raise RuntimeError(f"Only sampled {len(controls)} controls out of requested {n_controls}")
    return pd.DataFrame(controls)


def summarize(df, positives_initial, windows):
    rows = [
        {"metric": "n_rows", "value": len(df)},
        {"metric": "n_positive_initial_candidates", "value": positives_initial},
        {"metric": "n_positive_retained", "value": int(df["label"].eq(1).sum())},
        {"metric": "n_control_retained", "value": int(df["label"].eq(0).sum())},
        {"metric": "windows", "value": ",".join(map(str, windows))},
        {"metric": "junction_mask_bp", "value": df.attrs.get("junction_mask_bp", pd.NA)},
        {"metric": "n_reverse_oriented_rows", "value": int(df["reverse_oriented_contig"].fillna(0).astype(int).sum()) if "reverse_oriented_contig" in df else 0},
    ]
    positives = df.loc[df["label"].eq(1)]
    if not positives.empty:
        rows.extend([
            {"metric": "positive_median_absent_span", "value": positives["absent_span"].median()},
            {"metric": "positive_max_absent_span", "value": positives["absent_span"].max()},
            {"metric": "positive_n_span_le_200", "value": int(positives["absent_span"].le(200).sum())},
            {"metric": "positive_n_span_le_225", "value": int(positives["absent_span"].le(225).sum())},
        ])
    for label, group in df.groupby("label"):
        rows.append({"metric": f"label_{label}_n_genes", "value": group["gene"].nunique()})
        rows.append({"metric": f"label_{label}_n_contigs", "value": group["contig"].nunique()})
    return pd.DataFrame(rows)


def main():
    args = parse_args()
    windows = parse_windows(args.windows)
    max_window = max(windows)
    output = Path(args.output)
    summary = Path(args.summary)
    output.parent.mkdir(parents=True, exist_ok=True)
    summary.parent.mkdir(parents=True, exist_ok=True)

    print("Loading genotype matrix...")
    genotype_df = pd.read_csv(args.genotype_matrix, sep="\t", low_memory=False)
    genotype_df["start"] = genotype_df["start"].astype(int)
    genotype_df["end"] = genotype_df["end"].astype(int)

    print("Loading RCC1749 exons and FASTA...")
    exons_by_contig = read_exons(args.gtf)
    excluded_contigs = set(args.exclude_contig)
    for contig in excluded_contigs:
        exons_by_contig.pop(contig, None)
    fasta_sizes = read_fai(args.fa)
    fasta_index = SeqIO.index(args.fa, "fasta")
    args.reverse_contigs = read_reverse_oriented_contigs(
        args.mummer_coords,
        args.orientation_threshold,
    )
    if args.reverse_contigs:
        print(f"Orientation-normalizing {len(args.reverse_contigs)} reverse-oriented contigs.")

    initial_positive_ids = (
        genotype_df.loc[
            genotype_df["sample"].isin([args.reference_sample, args.sample]),
            ["ortholog_id", "sample", "presence"],
        ]
        .pivot_table(index="ortholog_id", columns="sample", values="presence", aggfunc="first")
    )
    n_initial = int(
        (
            initial_positive_ids[args.reference_sample].eq(1)
            & initial_positive_ids[args.sample].eq(2)
        ).sum()
    )

    print("Selecting positive pre-insertion sites...")
    positives = positive_rows(genotype_df, exons_by_contig, fasta_sizes, max_window, args)
    n_controls = len(positives) * args.controls_per_positive
    print(f"Retained {len(positives)} positives; sampling {n_controls} controls")

    print("Sampling exonic controls...")
    exclusion_index = build_exclusion_index(genotype_df, args.sample)
    controls = sample_control_rows(
        exons_by_contig,
        fasta_sizes,
        exclusion_index,
        n_controls,
        max_window,
        args.sample,
        args.random_state,
        reverse_contigs=args.reverse_contigs,
        excluded_contigs=excluded_contigs,
    )

    metadata_df = pd.concat([positives, controls], ignore_index=True)
    feature_rows = []
    print("Computing sequence features...")
    for row in metadata_df.to_dict("records"):
        feature_rows.append(
            add_window_features(
                row,
                fasta_index,
                windows,
                junction_mask_bp=args.junction_mask_bp,
            )
        )
    feature_df = pd.DataFrame(feature_rows)
    out_df = pd.concat([metadata_df, feature_df], axis=1)
    out_df["group_id"] = out_df["gene"].where(out_df["gene"].astype(str).ne(""), out_df["contig"])
    out_df.attrs["junction_mask_bp"] = args.junction_mask_bp

    out_df.to_csv(output, sep="\t", index=False)
    summarize(out_df, n_initial, windows).to_csv(summary, sep="\t", index=False)
    print(f"Wrote {output}")
    print(f"Wrote {summary}")


if __name__ == "__main__":
    main()
