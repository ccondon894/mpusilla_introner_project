import argparse
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


CODING_WINDOWS = [3, 9, 15, 30, 60]
VALID_BASES = set("ACGT")

GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "Stop", "TAG": "Stop",
    "TGT": "C", "TGC": "C", "TGA": "Stop", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

CODONS_BY_AA = defaultdict(list)
for codon, aa in GENETIC_CODE.items():
    CODONS_BY_AA[aa].append(codon)

TARGET_FAMILIES = {
    "pro": {"CCA", "CCC", "CCG", "CCT"},
    "gly": {"GGA", "GGC", "GGG", "GGT"},
    "ala": {"GCA", "GCC", "GCG", "GCT"},
    "arg": {"AGA", "AGG", "CGA", "CGC", "CGG", "CGT"},
    "leu": {"CTA", "CTC", "CTG", "CTT", "TTA", "TTG"},
    "ser": {"AGC", "AGT", "TCA", "TCC", "TCG", "TCT"},
    "gln": {"CAA", "CAG"},
}

DEGENERACY_LABELS = {
    1: "singleton",
    2: "twofold",
    3: "threefold",
    4: "fourfold",
    6: "sixfold",
}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--matrix",
        default="genotype_matrix_with_features.microc_tracks.tsv",
    )
    parser.add_argument(
        "--gtf",
        default="/scratch1/chris/mpusilla_introner_project/results/annotations/CCMP1545.gtf",
    )
    parser.add_argument(
        "--fa",
        default="/scratch1/alex/pusilla/data/ref/v3/MpusillaCCMP1545_228_v3.0.fa",
    )
    parser.add_argument(
        "--output",
        default="rf_results/genotype_matrix_with_features.microc_tracks.codon.tsv",
    )
    parser.add_argument(
        "--summary",
        default="rf_results/codon_feature_join_summary.tsv",
    )
    parser.add_argument(
        "--unmatched",
        default="rf_results/codon_feature_unmatched_loci.tsv",
    )
    parser.add_argument("--max-boundary-distance", type=int, default=100)
    parser.add_argument("--limit-rows", type=int, default=None)
    return parser.parse_args()


def parse_gtf_attrs(attr_text):
    return dict(re.findall(r'(\S+) "([^"]+)"', attr_text))


def fasta_key_for_contig(fasta_index, contig):
    if contig in fasta_index:
        return contig
    suffix = contig.split("#")[-1]
    if suffix in fasta_index:
        return suffix
    return None


def fetch_interval(fasta_index, contig, start_1based, end_1based, strand):
    fasta_key = fasta_key_for_contig(fasta_index, contig)
    if fasta_key is None or start_1based > end_1based:
        return ""
    seq = str(fasta_index[fasta_key].seq[start_1based - 1:end_1based].upper())
    if strand == "-":
        return str(Seq(seq).reverse_complement())
    return seq


def read_cds_by_transcript(gtf_path):
    cds_by_tx = defaultdict(list)
    with open(gtf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            attrs = parse_gtf_attrs(fields[8])
            transcript_id = attrs.get("transcript_id")
            gene_id = attrs.get("gene_id")
            if not transcript_id or not gene_id:
                continue
            cds_by_tx[transcript_id].append({
                "contig": fields[0],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6],
                "phase": fields[7],
                "gene": gene_id,
                "transcript_id": transcript_id,
            })
    return cds_by_tx


def build_cds_introns(cds_by_tx, fasta_index):
    introns = []
    for transcript_id, records in cds_by_tx.items():
        if len(records) < 2:
            continue

        strands = {record["strand"] for record in records}
        contigs = {record["contig"] for record in records}
        genes = {record["gene"] for record in records}
        if len(strands) != 1 or len(contigs) != 1 or len(genes) != 1:
            continue

        strand = records[0]["strand"]
        transcript_order = sorted(
            records,
            key=lambda record: record["start"],
            reverse=(strand == "-"),
        )

        coding_pieces = []
        offsets = {}
        cursor = 0
        for record in transcript_order:
            seq = fetch_interval(
                fasta_index,
                record["contig"],
                record["start"],
                record["end"],
                strand,
            )
            offsets[(record["start"], record["end"])] = (cursor, cursor + len(seq))
            coding_pieces.append(seq)
            cursor += len(seq)
        coding_seq = "".join(coding_pieces)

        for upstream, downstream in zip(transcript_order, transcript_order[1:]):
            lower, upper = sorted([upstream, downstream], key=lambda record: record["start"])
            intron_start = lower["end"] + 1
            intron_end = upper["start"] - 1
            if intron_start > intron_end:
                continue
            boundary_pos = offsets[(upstream["start"], upstream["end"])][1]
            introns.append({
                "gene": upstream["gene"],
                "transcript_id": transcript_id,
                "contig": upstream["contig"],
                "strand": strand,
                "start": intron_start,
                "end": intron_end,
                "boundary_pos": boundary_pos,
                "coding_seq": coding_seq,
            })
    return introns


def entropy(values):
    values = [value for value in values if value]
    if not values:
        return math.nan
    counts = Counter(values)
    total = sum(counts.values())
    return -sum((count / total) * math.log2(count / total) for count in counts.values())


def split_codons(seq):
    seq = seq.upper()
    codons = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        if len(codon) == 3 and set(codon) <= VALID_BASES:
            codons.append(codon)
    return codons


def codon_windows_around_boundary(coding_seq, boundary_pos, window_bp):
    codons = split_codons(coding_seq)
    if not codons:
        return [], [], []

    n_codons = max(window_bp // 3, 1)
    total_codons = len(codons)
    phase = boundary_pos % 3

    if phase == 0:
        upstream_end = boundary_pos // 3
    else:
        upstream_end = math.ceil(boundary_pos / 3)
    downstream_start = boundary_pos // 3

    upstream_end = min(max(upstream_end, 0), total_codons)
    downstream_start = min(max(downstream_start, 0), total_codons)

    upstream_start = max(0, upstream_end - n_codons)
    downstream_end = min(total_codons, downstream_start + n_codons)

    upstream_indices = list(range(upstream_start, upstream_end))
    downstream_indices = list(range(downstream_start, downstream_end))
    combined_indices = sorted(set(upstream_indices + downstream_indices))

    upstream_codons = [codons[i] for i in upstream_indices]
    downstream_codons = [codons[i] for i in downstream_indices]
    combined_codons = [codons[i] for i in combined_indices]
    return upstream_codons, downstream_codons, combined_codons


def aa_for_codon(codon):
    return GENETIC_CODE.get(codon) if codon else None


def degeneracy_for_codon(codon):
    aa = aa_for_codon(codon)
    if aa is None or aa == "Stop":
        return math.nan
    return len(CODONS_BY_AA[aa])


def codon_stats(codons):
    stats = {}
    if not codons:
        return {
            "n_codons": 0,
            "gc1": math.nan,
            "gc2": math.nan,
            "gc3": math.nan,
            "codon_entropy": math.nan,
            "aa_entropy": math.nan,
            "syn_family_entropy": math.nan,
            "singleton_fraction": math.nan,
            "twofold_fraction": math.nan,
            "fourfold_fraction": math.nan,
            "sixfold_fraction": math.nan,
            "gc_ending_synonymous_fraction": math.nan,
            **{f"{name}_family_fraction": math.nan for name in TARGET_FAMILIES},
        }

    stats["n_codons"] = len(codons)
    for pos in range(3):
        bases = [codon[pos] for codon in codons]
        stats[f"gc{pos + 1}"] = sum(base in {"G", "C"} for base in bases) / len(bases)

    aas = [aa_for_codon(codon) for codon in codons]
    syn_families = [aa for aa in aas if aa and aa != "Stop" and len(CODONS_BY_AA[aa]) > 1]
    stats["codon_entropy"] = entropy(codons)
    stats["aa_entropy"] = entropy(aas)
    stats["syn_family_entropy"] = entropy(syn_families)

    degeneracies = [degeneracy_for_codon(codon) for codon in codons]
    valid_degeneracies = [value for value in degeneracies if not math.isnan(value)]
    for size, label in DEGENERACY_LABELS.items():
        if size == 3:
            continue
        stats[f"{label}_fraction"] = (
            sum(value == size for value in valid_degeneracies) / len(valid_degeneracies)
            if valid_degeneracies else math.nan
        )

    synonymous_codons = [
        codon for codon in codons
        if aa_for_codon(codon) not in {None, "Stop"}
        and len(CODONS_BY_AA[aa_for_codon(codon)]) > 1
    ]
    stats["gc_ending_synonymous_fraction"] = (
        sum(codon[2] in {"G", "C"} for codon in synonymous_codons) / len(synonymous_codons)
        if synonymous_codons else math.nan
    )

    for name, family_codons in TARGET_FAMILIES.items():
        stats[f"{name}_family_fraction"] = (
            sum(codon in family_codons for codon in codons) / len(codons)
        )
    return stats


def add_prefixed_stats(features, prefix, stats):
    for name, value in stats.items():
        features[f"{prefix}_{name}"] = value


def boundary_codon(coding_seq, boundary_pos, side):
    if side == "upstream":
        if boundary_pos <= 0:
            return None
        codon_start = ((boundary_pos - 1) // 3) * 3
    else:
        if boundary_pos >= len(coding_seq):
            return None
        codon_start = (boundary_pos // 3) * 3
    codon = coding_seq[codon_start:codon_start + 3]
    if len(codon) != 3 or set(codon) - VALID_BASES:
        return None
    return codon


def add_boundary_codon_features(features, coding_seq, boundary_pos, side):
    codon = boundary_codon(coding_seq, boundary_pos, side)
    aa = aa_for_codon(codon)
    degeneracy = degeneracy_for_codon(codon)
    prefix = f"codon_boundary_{side}"
    features[f"{prefix}_triplet"] = codon
    features[f"{prefix}_aa"] = aa
    features[f"{prefix}_is_synonymous_family"] = (
        int(aa not in {None, "Stop"} and len(CODONS_BY_AA[aa]) > 1)
        if codon else math.nan
    )
    features[f"{prefix}_degeneracy"] = degeneracy
    features[f"{prefix}_gc3"] = (
        int(codon[2] in {"G", "C"}) if codon else math.nan
    )


def build_features_for_match(match, max_boundary_distance):
    features = {
        "codon_match_found": 1,
        "codon_transcript_id": match["transcript_id"],
        "codon_strand": match["strand"],
        "codon_cds_intron_start": match["start"],
        "codon_cds_intron_end": match["end"],
        "codon_start_delta": match["start_delta"],
        "codon_end_delta": match["end_delta"],
        "codon_abs_start_delta": abs(match["start_delta"]),
        "codon_abs_end_delta": abs(match["end_delta"]),
        "codon_boundary_distance_sum": match["distance_sum"],
        "codon_boundary_match_type": (
            "exact" if match["distance_sum"] == 0 else
            "within_threshold" if match["distance_sum"] <= max_boundary_distance else
            "nearest_outside_threshold"
        ),
        "codon_has_left_cds_boundary": int(match["distance_sum"] <= max_boundary_distance),
        "codon_has_right_cds_boundary": int(match["distance_sum"] <= max_boundary_distance),
        "codon_has_both_cds_boundaries": int(match["distance_sum"] <= max_boundary_distance),
        "codon_intron_phase": match["boundary_pos"] % 3,
    }

    coding_seq = match["coding_seq"]
    boundary_pos = match["boundary_pos"]
    add_boundary_codon_features(features, coding_seq, boundary_pos, "upstream")
    add_boundary_codon_features(features, coding_seq, boundary_pos, "downstream")

    for window in CODING_WINDOWS:
        upstream_codons, downstream_codons, combined_codons = (
            codon_windows_around_boundary(coding_seq, boundary_pos, window)
        )

        upstream_stats = codon_stats(upstream_codons)
        downstream_stats = codon_stats(downstream_codons)
        combined_stats = codon_stats(combined_codons)

        add_prefixed_stats(features, f"codon_upstream_{window}", upstream_stats)
        add_prefixed_stats(features, f"codon_downstream_{window}", downstream_stats)
        add_prefixed_stats(features, f"codon_combined_{window}", combined_stats)

        for name, upstream_value in upstream_stats.items():
            downstream_value = downstream_stats[name]
            if name == "n_codons":
                features[f"codon_delta_{window}_{name}"] = upstream_value - downstream_value
            elif pd.notna(upstream_value) and pd.notna(downstream_value):
                features[f"codon_delta_{window}_{name}"] = upstream_value - downstream_value
            else:
                features[f"codon_delta_{window}_{name}"] = math.nan

    return features


def empty_features(reason):
    return {
        "codon_match_found": 0,
        "codon_boundary_match_type": reason,
        "codon_has_left_cds_boundary": 0,
        "codon_has_right_cds_boundary": 0,
        "codon_has_both_cds_boundaries": 0,
    }


def best_match(row, introns_by_gene_contig):
    gene = row.get("gene")
    contig = row.get("contig")
    if pd.isna(gene) or pd.isna(contig):
        return None, "missing_gene_or_contig"

    candidates = introns_by_gene_contig.get((str(gene), str(contig)), [])
    if not candidates:
        return None, "no_same_gene_cds_intron"

    row_start = int(row.get("feature_start", row["start"]))
    row_end = int(row.get("feature_end", row["end"]))
    scored = []
    for candidate in candidates:
        start_delta = row_start - candidate["start"]
        end_delta = row_end - candidate["end"]
        distance_sum = abs(start_delta) + abs(end_delta)
        scored.append((distance_sum, abs(start_delta), abs(end_delta), candidate, start_delta, end_delta))
    scored.sort(key=lambda item: (item[0], item[1], item[2], item[3]["transcript_id"]))
    best = scored[0]
    tied = [item for item in scored if item[:3] == best[:3]]
    if len(tied) > 1:
        return None, "ambiguous_transcript"

    match = dict(best[3])
    match["start_delta"] = best[4]
    match["end_delta"] = best[5]
    match["distance_sum"] = best[0]
    return match, None


def summarize_join(df):
    rows = []
    rows.append({"metric": "n_rows", "value": len(df)})
    rows.append({"metric": "n_matched", "value": int(df["codon_match_found"].sum())})
    rows.append({
        "metric": "n_within_100bp",
        "value": int((df["codon_boundary_distance_sum"] <= 100).sum()),
    })
    rows.append({
        "metric": "n_exact",
        "value": int((df["codon_boundary_distance_sum"] == 0).sum()),
    })
    for label, count in df["label"].value_counts(dropna=False).items():
        rows.append({"metric": f"n_label_{label}", "value": int(count)})
    for match_type, count in df["codon_boundary_match_type"].value_counts(dropna=False).items():
        rows.append({"metric": f"match_type_{match_type}", "value": int(count)})
    for label, group in df.groupby("label", dropna=False):
        rows.append({
            "metric": f"label_{label}_n_matched",
            "value": int(group["codon_match_found"].sum()),
        })
        rows.append({
            "metric": f"label_{label}_median_distance_sum",
            "value": group["codon_boundary_distance_sum"].median(),
        })
    return pd.DataFrame(rows)


def main():
    args = parse_args()
    output = Path(args.output)
    summary = Path(args.summary)
    unmatched = Path(args.unmatched)
    output.parent.mkdir(parents=True, exist_ok=True)
    summary.parent.mkdir(parents=True, exist_ok=True)
    unmatched.parent.mkdir(parents=True, exist_ok=True)

    print("Loading matrix...")
    matrix_df = pd.read_csv(args.matrix, sep="\t")
    if args.limit_rows is not None:
        matrix_df = matrix_df.head(args.limit_rows).copy()

    print("Indexing fasta...")
    fasta_index = SeqIO.index(args.fa, "fasta")

    print("Parsing GTF CDS features...")
    cds_by_tx = read_cds_by_transcript(args.gtf)
    cds_introns = build_cds_introns(cds_by_tx, fasta_index)
    introns_by_gene_contig = defaultdict(list)
    for intron in cds_introns:
        introns_by_gene_contig[(intron["gene"], intron["contig"])].append(intron)

    print("Adding codon features...")
    feature_rows = []
    for row in matrix_df.to_dict("records"):
        match, reason = best_match(row, introns_by_gene_contig)
        if match is None:
            feature_rows.append(empty_features(reason))
        else:
            feature_rows.append(build_features_for_match(match, args.max_boundary_distance))

    codon_df = pd.DataFrame(feature_rows)
    out_df = pd.concat([matrix_df.reset_index(drop=True), codon_df], axis=1)
    out_df.to_csv(output, sep="\t", index=False)

    summarize_join(out_df).to_csv(summary, sep="\t", index=False)

    unmatched_cols = [
        "sequence_id", "gene", "contig", "start", "end",
        "feature_start", "feature_end", "label",
        "codon_boundary_match_type", "codon_boundary_distance_sum",
        "codon_start_delta", "codon_end_delta",
    ]
    unmatched_cols = [col for col in unmatched_cols if col in out_df.columns]
    out_df.loc[out_df["codon_match_found"] == 0, unmatched_cols].to_csv(
        unmatched,
        sep="\t",
        index=False,
    )

    print(f"Wrote {output}")
    print(f"Wrote {summary}")
    print(f"Wrote {unmatched}")


if __name__ == "__main__":
    main()
