#!/usr/bin/env python3
"""Calculate gene-level GC content for expression GLM covariates."""

import argparse
import re
from collections import defaultdict

import pandas as pd


def parse_attributes(attr_string):
    attrs = {}
    for item in attr_string.split(";"):
        item = item.strip()
        if not item:
            continue
        match = re.match(r'([^=\s]+)[=\s]+"?([^";]+)"?', item)
        if match:
            attrs[match.group(1)] = match.group(2)
    return attrs


def load_fasta(path):
    seqs = {}
    name = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def reverse_complement(seq):
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def longest_exon_isoforms(gtf_path):
    transcripts = defaultdict(lambda: defaultdict(list))
    with open(gtf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            attrs = parse_attributes(fields[8])
            gene_id = attrs.get("gene_id")
            transcript_id = attrs.get("transcript_id")
            if not gene_id or not transcript_id:
                continue
            transcripts[gene_id][transcript_id].append({
                "chrom": fields[0],
                "start": int(fields[3]) - 1,
                "end": int(fields[4]),
                "strand": fields[6],
            })

    longest = {}
    for gene_id, by_transcript in transcripts.items():
        best_id = None
        best_exons = None
        best_len = -1
        for transcript_id, exons in by_transcript.items():
            total_len = sum(exon["end"] - exon["start"] for exon in exons)
            if total_len > best_len:
                best_id = transcript_id
                best_exons = exons
                best_len = total_len
        if best_id is not None:
            longest[gene_id] = sorted(best_exons, key=lambda exon: exon["start"])
    return longest


def gc_fraction(seq):
    seq = seq.upper()
    called = sum(base in "ACGT" for base in seq)
    if called == 0:
        return None
    return (seq.count("G") + seq.count("C")) / called


def calculate_for_sample(sample, gtf_path, fasta_path):
    fasta = load_fasta(fasta_path)
    rows = []
    for gene_id, exons in longest_exon_isoforms(gtf_path).items():
        pieces = []
        missing = False
        for exon in exons:
            chrom_seq = fasta.get(exon["chrom"])
            if chrom_seq is None:
                missing = True
                break
            pieces.append(chrom_seq[exon["start"]:exon["end"]])
        if missing or not pieces:
            continue
        seq = "".join(pieces)
        if exons[0]["strand"] == "-":
            seq = reverse_complement(seq)
        gc = gc_fraction(seq)
        if gc is not None:
            rows.append({"gene_id": gene_id, "strain": sample, "GC_content": gc})
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", nargs="+", required=True)
    parser.add_argument("--gtfs", nargs="+", required=True)
    parser.add_argument("--assemblies", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    if not (len(args.samples) == len(args.gtfs) == len(args.assemblies)):
        raise ValueError("--samples, --gtfs, and --assemblies must have equal lengths")

    rows = []
    for sample, gtf_path, fasta_path in zip(args.samples, args.gtfs, args.assemblies):
        sample_rows = calculate_for_sample(sample, gtf_path, fasta_path)
        print(f"{sample}: calculated GC content for {len(sample_rows)} genes")
        rows.extend(sample_rows)

    if not rows:
        raise ValueError("No GC content values were calculated")

    pd.DataFrame(rows).to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
