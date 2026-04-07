#!/usr/bin/env python3
"""
Diagnostic: check whether introns polarized as "gains" are actually present
in outgroup assemblies (i.e., false negatives from annotation/alignment).

Reads the genotype matrix, identifies introns where both outgroup samples
are ABSENT (polarized as gains), extracts reference sequences, and searches
them against outgroup assemblies using minimap2 with sensitive settings.
"""

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pysam

PRESENT = 1
ABSENT = 2
MISSING = 3


def main():
    parser = argparse.ArgumentParser(description="Diagnose gain introns against outgroup assemblies")
    parser.add_argument("--matrix", required=True, help="Intron genotype matrix TSV")
    parser.add_argument("--assemblies_dir", required=True, help="Directory with {sample}.vg_paths.fa")
    parser.add_argument("--reference", default="CCMP1545", help="Reference sample")
    parser.add_argument("--group2", required=True, help="Comma-separated outgroup samples")
    parser.add_argument("--output", required=True, help="Output diagnostic report")
    args = parser.parse_args()

    group2 = args.group2.split(",")

    # Parse matrix to find gain introns
    # A "gain" = both outgroup ABSENT, at least one Group 1 sample PRESENT
    gains = []
    with open(args.matrix) as f:
        header = f.readline().rstrip("\n").split("\t")
        sample_cols = header[5:]  # after gene_id, contig, ref_start, ref_end, intron_index

        for line in f:
            fields = line.rstrip("\n").split("\t")
            gene_id = fields[0]
            contig = fields[1]
            ref_start = int(fields[2])
            ref_end = int(fields[3])
            intron_idx = fields[4]

            genotypes = {s: int(v) for s, v in zip(sample_cols, fields[5:])}

            # Check outgroup: all ABSENT (not MISSING)
            g2_states = [genotypes.get(s) for s in group2]
            if any(s != ABSENT for s in g2_states):
                continue

            # Only include polymorphic gains (matches unfolded AFS criteria)
            # Exclude sites with MISSING in Group 1
            g1_samples = [s for s in sample_cols if s not in group2]
            g1_states = [genotypes.get(s) for s in g1_samples]
            if MISSING in g1_states:
                continue
            # Derived = present (since ancestral = absent)
            derived_count = g1_states.count(PRESENT)
            # Skip fixed (DAC=0 or DAC=n)
            if derived_count == 0 or derived_count == len(g1_samples):
                continue

            gains.append({
                "gene_id": gene_id,
                "contig": contig,
                "ref_start": ref_start,
                "ref_end": ref_end,
                "intron_idx": intron_idx,
                "genotypes": genotypes,
            })

    print(f"Found {len(gains)} introns polarized as gains")

    # Extract reference sequences for gain introns
    ref_fasta_path = os.path.join(args.assemblies_dir, f"{args.reference}.vg_paths.fa")
    ref_fa = pysam.FastaFile(ref_fasta_path)

    tmpdir = tempfile.mkdtemp(prefix="diagnose_gains_")
    query_fasta = os.path.join(tmpdir, "gain_introns.fa")

    with open(query_fasta, "w") as f:
        for g in gains:
            seq = ref_fa.fetch(g["contig"], g["ref_start"] - 1, g["ref_end"])
            g["sequence"] = seq
            g["length"] = len(seq)
            name = f"{g['gene_id']}__{g['intron_idx']}__{g['contig']}:{g['ref_start']}-{g['ref_end']}"
            f.write(f">{name}\n{seq}\n")

    ref_fa.close()

    # Search each outgroup assembly
    outgroup_hits = {s: {} for s in group2}

    for sample in group2:
        assembly = os.path.join(args.assemblies_dir, f"{sample}.vg_paths.fa")
        paf_path = os.path.join(tmpdir, f"{sample}.paf")

        # Use very sensitive minimap2 settings for short sequences
        cmd = [
            "minimap2", "-k", "8", "-w", "3", "--secondary=no", "-c",
            "-N", "5",  # report up to 5 hits
            assembly, query_fasta,
        ]
        with open(paf_path, "w") as paf_out:
            subprocess.run(cmd, stdout=paf_out, stderr=subprocess.DEVNULL, check=True)

        # Parse PAF
        with open(paf_path) as f:
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 12:
                    continue
                query_name = cols[0]
                query_len = int(cols[1])
                target_name = cols[5]
                nmatch = int(cols[9])
                block_len = int(cols[10])
                mapq = int(cols[11])

                identity = nmatch / block_len if block_len > 0 else 0
                query_cov = nmatch / query_len if query_len > 0 else 0

                # Keep best hit per query
                key = query_name
                if key not in outgroup_hits[sample] or identity > outgroup_hits[sample][key]["identity"]:
                    outgroup_hits[sample][key] = {
                        "target": target_name,
                        "identity": identity,
                        "query_cov": query_cov,
                        "nmatch": nmatch,
                        "block_len": block_len,
                        "mapq": mapq,
                    }

        os.remove(paf_path)

    os.remove(query_fasta)
    os.rmdir(tmpdir)

    # Write report
    with open(args.output, "w") as f:
        f.write("=== Gain Intron Diagnostic Report ===\n\n")

        # Summary counts
        n_found_any = 0
        n_found_all = 0
        n_found_none = 0

        for g in gains:
            name = f"{g['gene_id']}__{g['intron_idx']}__{g['contig']}:{g['ref_start']}-{g['ref_end']}"
            found_in = []
            for s in group2:
                hit = outgroup_hits[s].get(name)
                if hit and hit["identity"] > 0.5:
                    found_in.append(s)
            if len(found_in) == len(group2):
                n_found_all += 1
            elif len(found_in) > 0:
                n_found_any += 1
            else:
                n_found_none += 1

        f.write(f"Total gain introns: {len(gains)}\n")
        f.write(f"Found in ALL outgroup assemblies (false negatives): {n_found_all}\n")
        f.write(f"Found in SOME outgroup assemblies: {n_found_any}\n")
        f.write(f"Found in NO outgroup assemblies (possibly real gains): {n_found_none}\n\n")

        # Per-intron details
        f.write(f"{'gene_id':<20} {'idx':>3} {'len':>5} {'contig':<25} {'start':>8} {'end':>8}")
        for s in group2:
            f.write(f"  {s + '_id':>12} {s + '_cov':>10}")
        f.write(f"  {'verdict':<15}\n")
        f.write("-" * (90 + 24 * len(group2)) + "\n")

        for g in gains:
            name = f"{g['gene_id']}__{g['intron_idx']}__{g['contig']}:{g['ref_start']}-{g['ref_end']}"
            f.write(f"{g['gene_id']:<20} {g['intron_idx']:>3} {g['length']:>5} {g['contig']:<25} {g['ref_start']:>8} {g['ref_end']:>8}")

            found_in = []
            for s in group2:
                hit = outgroup_hits[s].get(name)
                if hit:
                    f.write(f"  {hit['identity']:>12.3f} {hit['query_cov']:>10.3f}")
                    if hit["identity"] > 0.5:
                        found_in.append(s)
                else:
                    f.write(f"  {'no_hit':>12} {'':>10}")

            if len(found_in) == len(group2):
                verdict = "FALSE_NEG"
            elif len(found_in) > 0:
                verdict = "PARTIAL"
            else:
                verdict = "REAL_GAIN?"
            f.write(f"  {verdict:<15}\n")

    print(f"\nReport written to: {args.output}")
    print(f"  Found in all outgroups (false negatives): {n_found_all}")
    print(f"  Found in some outgroups: {n_found_any}")
    print(f"  Found in no outgroups (possibly real): {n_found_none}")


if __name__ == "__main__":
    main()
