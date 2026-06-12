#!/usr/bin/env python3
"""
Non-introner intron presence/absence polymorphism and allele frequency spectrum.

Extracts introns from GTFs across all samples, filters out introner introns,
uses minimap2 sequence alignment to identify orthologous introns across samples,
builds a presence/absence genotype matrix, and computes both folded and unfolded
allele frequency spectra.

v2: Uses minimap2 alignment instead of cumulative CDS position matching to avoid
false positives from miniprot exon boundary variation across assemblies.
"""

import argparse
import os
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pysam


def save_figure_with_png(output_path):
    """Save the current matplotlib figure and emit a sibling PNG for PDF outputs."""
    output_path = Path(output_path)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    if output_path.suffix.lower() == ".pdf":
        plt.savefig(output_path.with_suffix(".png"), dpi=300, bbox_inches="tight")


# ── Stage 1: Extract introns from GTF ─────────────────────────────────────────

def parse_gtf_attribute(attr_str, key="transcript_id"):
    """Extract a value from a GTF attribute string."""
    for field in attr_str.strip().split(";"):
        field = field.strip()
        if field.startswith(key):
            return field.split('"')[1]
    return None


def extract_introns_from_gtf(gtf_path):
    """
    Parse a GTF and extract introns for each gene.

    Returns dict: gene_id -> list of (contig, start, end, intron_index)
    where intron_index is 0-based ordinal in coding order
    (plus strand: genomic order; minus strand: reversed).
    """
    gene_exons = defaultdict(list)
    gene_strand = {}
    gene_contig = {}

    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != "exon":
                continue

            contig = parts[0]
            start = int(parts[3])  # 1-based inclusive
            end = int(parts[4])    # 1-based inclusive
            strand = parts[6]
            gene_id = parse_gtf_attribute(parts[8])
            if gene_id is None:
                continue

            gene_exons[gene_id].append((start, end))
            gene_strand[gene_id] = strand
            gene_contig[gene_id] = contig

    introns = {}
    for gene_id, exons in gene_exons.items():
        if len(exons) < 2:
            continue

        strand = gene_strand[gene_id]
        contig = gene_contig[gene_id]

        # Sort exons by genomic position
        exons_sorted = sorted(exons, key=lambda x: x[0])

        # Compute intron gaps in genomic order
        genomic_introns = []
        for i in range(len(exons_sorted) - 1):
            intron_start = exons_sorted[i][1] + 1
            intron_end = exons_sorted[i + 1][0] - 1
            genomic_introns.append((intron_start, intron_end))

        # Assign coding-order indices
        if strand == "-":
            # Minus strand: coding order is reversed genomic order
            result = []
            for idx, (istart, iend) in enumerate(reversed(genomic_introns)):
                result.append((contig, istart, iend, idx))
        else:
            result = []
            for idx, (istart, iend) in enumerate(genomic_introns):
                result.append((contig, istart, iend, idx))

        introns[gene_id] = result

    return introns


# ── Stage 2: Extract intron sequences from assemblies ──────────────────────────

def extract_intron_sequences(sample_introns, assemblies_dir, sample):
    """
    Extract intron DNA sequences from assembly FASTA using pysam.

    Returns dict: gene_id -> list of (intron_index, sequence)
    """
    fasta_path = os.path.join(assemblies_dir, f"{sample}.vg_paths.fa")
    fa = pysam.FastaFile(fasta_path)

    sequences = {}
    for gene_id, intron_list in sample_introns.items():
        gene_seqs = []
        for contig, start, end, intron_idx in intron_list:
            # pysam.fetch uses 0-based half-open; GTF coords are 1-based inclusive
            try:
                seq = fa.fetch(contig, start - 1, end)
            except (KeyError, ValueError):
                continue
            if seq:
                gene_seqs.append((intron_idx, seq))
        if gene_seqs:
            sequences[gene_id] = gene_seqs

    fa.close()
    return sequences


# ── Stage 3: Filter introner introns from reference ────────────────────────────

def load_introner_loci(bed_path):
    """
    Load introner loci BED file and return set of (contig, body_start, body_end).
    The BED coordinates include 100bp flanking on each side, so the actual
    introner body is (start + 100, end - 100).
    """
    loci = []
    with open(bed_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            contig = parts[0]
            bed_start = int(parts[1])  # 0-based
            bed_end = int(parts[2])    # 0-based exclusive
            # Strip 100bp flanking on each side to get introner body
            body_start = bed_start + 100
            body_end = bed_end - 100
            loci.append((contig, body_start, body_end))
    return loci


def overlaps(intron_start, intron_end, locus_start, locus_end):
    """Check if two intervals overlap. intron coords are 1-based inclusive,
    locus coords are 0-based half-open from BED."""
    # Convert intron to 0-based half-open for comparison
    return intron_start - 1 < locus_end and locus_start < intron_end


def filter_introner_introns(ref_introns, introner_loci):
    """
    Remove reference introns that overlap with introner loci.
    Returns filtered dict of gene_id -> list of (contig, start, end, intron_index).
    """
    loci_by_contig = defaultdict(list)
    for contig, body_start, body_end in introner_loci:
        loci_by_contig[contig].append((body_start, body_end))

    filtered = {}
    n_removed = 0
    n_kept = 0

    for gene_id, intron_list in ref_introns.items():
        kept = []
        for contig, istart, iend, intron_idx in intron_list:
            is_introner = False
            for lstart, lend in loci_by_contig.get(contig, []):
                if overlaps(istart, iend, lstart, lend):
                    is_introner = True
                    break
            if is_introner:
                n_removed += 1
            else:
                kept.append((contig, istart, iend, intron_idx))
                n_kept += 1
        if kept:
            filtered[gene_id] = kept

    return filtered, n_removed, n_kept


# ── Stage 3b: Load coverage calls ─────────────────────────────────────────────

def load_coverage_calls(coverage_dir, samples):
    """
    Load coverage calling results for all samples.
    Returns dict: (sample, intron_id) -> new_call string (e.g. "1", "2", "3")
    """
    calls = {}
    for sample in samples:
        path = os.path.join(coverage_dir, f"{sample}.intron_coverage_calls.tsv")
        if not os.path.exists(path):
            print(f"  Warning: coverage calls not found for {sample}: {path}")
            continue
        n = 0
        with open(path) as f:
            header = f.readline()
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 10:
                    continue
                intron_id = parts[3]  # intron_id column
                new_call = parts[9]   # new_call column
                calls[(sample, intron_id)] = new_call
                n += 1
        print(f"  Loaded {n} coverage calls for {sample}")
    return calls


# ── Stage 4: minimap2 ortholog matching ────────────────────────────────────────

def run_minimap2_matching(ref_nonintroner, ref_sequences, all_intron_seqs,
                          all_sample_introns, group1, group2, reference,
                          coverage_calls=None):
    """
    Use minimap2 to match each sample's introns against reference non-introner introns.
    When coverage_calls is provided, downgrades PRESENT to ABSENT if coverage says absent/missing.

    Returns dict: (gene_id, intron_idx) -> {sample: PRESENT/ABSENT/MISSING}
    """
    all_samples = group1 + group2

    # Build set of reference non-introner intron keys
    ref_intron_keys = set()
    for gene_id, intron_list in ref_nonintroner.items():
        for contig, istart, iend, intron_idx in intron_list:
            ref_intron_keys.add((gene_id, intron_idx))

    # Write reference non-introner intron sequences to temp FASTA
    tmpdir = tempfile.mkdtemp(prefix="minimap2_introns_")
    ref_fasta = os.path.join(tmpdir, "ref_introns.fa")

    n_ref_written = 0
    with open(ref_fasta, "w") as f:
        for gene_id, intron_idx in sorted(ref_intron_keys):
            # Find the sequence from ref_sequences
            if gene_id not in ref_sequences:
                continue
            for idx, seq in ref_sequences[gene_id]:
                if idx == intron_idx:
                    f.write(f">{gene_id}__intron_{intron_idx}\n{seq}\n")
                    n_ref_written += 1
                    break

    print(f"  Wrote {n_ref_written} reference non-introner intron sequences")

    # Initialize results: reference has all its non-introner introns as PRESENT
    results = {key: {} for key in ref_intron_keys}
    for key in ref_intron_keys:
        results[key][reference] = PRESENT

    # For each non-reference sample, run minimap2
    for sample in all_samples:
        if sample == reference:
            continue

        sample_introns = all_sample_introns.get(sample, {})
        sample_seqs = all_intron_seqs.get(sample, {})

        # Genes present in this sample (from GTF)
        sample_genes = set(sample_introns.keys())

        # Write sample intron sequences (only for genes in reference)
        sample_fasta = os.path.join(tmpdir, f"{sample}_introns.fa")
        n_sample_written = 0
        ref_genes = set(g for g, _ in ref_intron_keys)

        with open(sample_fasta, "w") as f:
            for gene_id in ref_genes:
                if gene_id not in sample_seqs:
                    continue
                for intron_idx, seq in sample_seqs[gene_id]:
                    f.write(f">{gene_id}__intron_{intron_idx}\n{seq}\n")
                    n_sample_written += 1

        if n_sample_written == 0:
            # No sequences to align — mark all as MISSING or ABSENT
            for key in ref_intron_keys:
                gene_id, _ = key
                if gene_id not in sample_genes:
                    results[key][sample] = MISSING
                else:
                    results[key][sample] = ABSENT
            os.remove(sample_fasta)
            continue

        # Run minimap2
        paf_path = os.path.join(tmpdir, f"{sample}.paf")
        cmd = [
            "minimap2", "-k", "10", "-w", "5", "--secondary=no", "-c",
            ref_fasta, sample_fasta,
        ]
        with open(paf_path, "w") as paf_out:
            subprocess.run(cmd, stdout=paf_out, stderr=subprocess.DEVNULL, check=True)

        # Parse PAF: find matching reference introns
        matched_ref_introns = set()
        with open(paf_path) as f:
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 12:
                    continue

                query_name = cols[0]   # sample intron
                target_name = cols[5]  # reference intron
                nmatch = int(cols[9])  # number of matching bases
                block_len = int(cols[10])  # alignment block length

                # Extract gene_ids from headers
                q_gene = query_name.rsplit("__intron_", 1)[0]
                t_gene = target_name.rsplit("__intron_", 1)[0]

                # Only accept within-gene matches
                if q_gene != t_gene:
                    continue

                # Check identity
                identity = nmatch / block_len if block_len > 0 else 0
                if identity > 0.5:
                    # Parse reference intron index from target name
                    t_idx = int(target_name.rsplit("__intron_", 1)[1])
                    matched_ref_introns.add((t_gene, t_idx))

        # Assign genotypes for this sample
        n_cov_downgraded = 0
        for key in ref_intron_keys:
            gene_id, intron_idx = key
            if gene_id not in sample_genes:
                results[key][sample] = MISSING
            elif key in matched_ref_introns:
                # Coverage validation: only downgrade if coverage confidently says absent (2)
                # If coverage is missing (3) or unavailable, trust minimap2
                if coverage_calls is not None:
                    intron_id = f"{gene_id}__intron_{intron_idx}"
                    cov_call = coverage_calls.get((sample, intron_id))
                    if cov_call == "2":
                        results[key][sample] = ABSENT
                        n_cov_downgraded += 1
                    else:
                        results[key][sample] = PRESENT
                else:
                    results[key][sample] = PRESENT
            else:
                results[key][sample] = ABSENT

        # Clean up sample files
        os.remove(sample_fasta)
        os.remove(paf_path)
        cov_msg = f", {n_cov_downgraded} coverage-downgraded" if coverage_calls is not None else ""
        print(f"  {sample}: {n_sample_written} introns aligned, "
              f"{len(matched_ref_introns)} ref introns matched{cov_msg}")

    # Clean up
    os.remove(ref_fasta)
    os.rmdir(tmpdir)

    return results


# ── Stage 5: Build matrix, compute AFS, plot ──────────────────────────────────

PRESENT = 1
ABSENT = 2
MISSING = 3


def build_genotype_matrix(ref_nonintroner, genotype_results, group1, group2):
    """
    Build presence/absence matrix from minimap2 matching results.
    Returns list of dicts with intron info and genotype per sample.
    """
    all_samples = group1 + group2
    rows = []

    for gene_id, intron_list in sorted(ref_nonintroner.items()):
        for contig, istart, iend, intron_idx in intron_list:
            key = (gene_id, intron_idx)
            if key not in genotype_results:
                continue

            row = {
                "gene_id": gene_id,
                "contig": contig,
                "ref_start": istart,
                "ref_end": iend,
                "intron_index": intron_idx,
            }
            for sample in all_samples:
                row[sample] = genotype_results[key].get(sample, MISSING)

            rows.append(row)

    return rows


def compute_folded_afs(matrix, group1):
    """
    Compute folded AFS for Group 1 samples.
    Only uses introns with complete data (no MISSING) in Group 1.
    Minor allele frequency = min(present_count, n - present_count).
    Range: 1 to n//2 (since fixed sites are excluded from the spectrum).

    Also splits into minor=present vs minor=absent categories.
    """
    n = len(group1)
    afs = defaultdict(int)
    afs_minor_present = defaultdict(int)  # minor allele = intron present (rare intron)
    afs_minor_absent = defaultdict(int)   # minor allele = intron absent (rare loss)

    for row in matrix:
        genotypes = [row[s] for s in group1]
        if MISSING in genotypes:
            continue
        present_count = genotypes.count(PRESENT)
        # Skip fixed sites (all present or all absent)
        if present_count == 0 or present_count == n:
            continue
        maf = min(present_count, n - present_count)
        afs[maf] += 1
        if present_count <= n - present_count:
            afs_minor_present[maf] += 1
        else:
            afs_minor_absent[maf] += 1

    return dict(afs), dict(afs_minor_present), dict(afs_minor_absent), n


def compute_unfolded_afs(matrix, group1, group2):
    """
    Compute unfolded AFS using Group 2 as outgroup.
    Ancestral state inferred when both Group 2 samples agree.
    Derived allele count in Group 1 ranges 1 to n-1.
    Skip introns where Group 2 is polymorphic, missing, or Group 1 has missing.

    Returns separate loss and gain AFS dicts in addition to combined.
    """
    n = len(group1)
    afs = defaultdict(int)
    afs_loss = defaultdict(int)   # derived = intron loss
    afs_gain = defaultdict(int)   # derived = intron gain
    n_skipped_outgroup = 0

    for row in matrix:
        g1_genotypes = [row[s] for s in group1]
        g2_genotypes = [row[s] for s in group2]

        # Skip if any Group 1 missing
        if MISSING in g1_genotypes:
            continue
        # Skip if any Group 2 missing
        if MISSING in g2_genotypes:
            n_skipped_outgroup += 1
            continue
        # Skip if Group 2 is polymorphic (not all agree)
        if len(set(g2_genotypes)) > 1:
            n_skipped_outgroup += 1
            continue

        ancestral = g2_genotypes[0]  # all agree
        present_count = g1_genotypes.count(PRESENT)

        # Derived allele count
        if ancestral == PRESENT:
            # Ancestral = present → derived = absent (intron loss)
            derived_count = n - present_count
            target = afs_loss
        else:
            # Ancestral = absent → derived = present (intron gain)
            derived_count = present_count
            target = afs_gain

        # Skip fixed ancestral (derived_count=0) and fixed derived (derived_count=n)
        if derived_count == 0 or derived_count == n:
            continue

        afs[derived_count] += 1
        target[derived_count] += 1

    return dict(afs), dict(afs_loss), dict(afs_gain), n, n_skipped_outgroup


def _plot_unfolded_panel(ax, afs_loss, afs_gain, n, title):
    """Plot a single unfolded AFS panel (stacked loss/gain)."""
    x = list(range(1, n))
    y_loss = [afs_loss.get(i, 0) for i in x]
    y_gain = [afs_gain.get(i, 0) for i in x]
    total_loss = sum(y_loss)
    total_gain = sum(y_gain)

    ax.bar(x, y_loss, color="#E24A33", edgecolor="black", linewidth=0.5,
           label=f"Intron loss ({total_loss})")
    ax.bar(x, y_gain, bottom=y_loss, color="#6ACC65", edgecolor="black",
           linewidth=0.5, label=f"Intron gain ({total_gain})")
    ax.set_xlabel(f"Derived allele count (Group 1, n={n})")
    ax.set_ylabel("Number of intron loci")
    ax.set_title(title)
    ax.set_xticks(x)
    ax.legend(loc="upper right", fontsize=9)


def plot_afs(folded_afs_minor_present, folded_afs_minor_absent, folded_n,
             unfolded_afs_loss, unfolded_afs_gain, unfolded_n,
             unfolded_1749_afs_loss, unfolded_1749_afs_gain, unfolded_1749_n,
             output_path):
    """Create three-panel figure: folded, unfolded (all Group 2), unfolded (RCC1749 only)."""
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(7, 11))

    # --- Folded AFS (stacked minor=present / minor=absent) ---
    max_maf = folded_n // 2
    x_folded = list(range(1, max_maf + 1))
    y_minor_present = [folded_afs_minor_present.get(i, 0) for i in x_folded]
    y_minor_absent = [folded_afs_minor_absent.get(i, 0) for i in x_folded]
    total_minor_present = sum(y_minor_present)
    total_minor_absent = sum(y_minor_absent)

    ax1.bar(x_folded, y_minor_absent, color="#E24A33", edgecolor="black", linewidth=0.5,
            label=f"Minor allele = absent ({total_minor_absent})")
    ax1.bar(x_folded, y_minor_present, bottom=y_minor_absent, color="#6ACC65",
            edgecolor="black", linewidth=0.5,
            label=f"Minor allele = present ({total_minor_present})")
    ax1.set_xlabel("Minor allele count (Group 1, n=11)")
    ax1.set_ylabel("Number of intron loci")
    ax1.set_title("Folded AFS — Non-introner introns")
    ax1.set_xticks(x_folded)
    ax1.legend(loc="upper right", fontsize=9)

    # --- Unfolded AFS: all Group 2 as outgroup ---
    _plot_unfolded_panel(ax2, unfolded_afs_loss, unfolded_afs_gain, unfolded_n,
                         "Unfolded AFS — Outgroup: Group 2 (RCC1749 + RCC3052)")

    # --- Unfolded AFS: RCC1749 only as outgroup ---
    _plot_unfolded_panel(ax3, unfolded_1749_afs_loss, unfolded_1749_afs_gain, unfolded_1749_n,
                         "Unfolded AFS — Outgroup: RCC1749 only")

    plt.tight_layout()
    save_figure_with_png(output_path)
    plt.close()


# ── Output functions ───────────────────────────────────────────────────────────

def write_matrix(matrix, group1, group2, output_path):
    """Write genotype matrix to TSV."""
    samples = group1 + group2
    header = ["gene_id", "contig", "ref_start", "ref_end", "intron_index"] + samples

    with open(output_path, "w") as f:
        f.write("\t".join(header) + "\n")
        for row in matrix:
            fields = [
                row["gene_id"],
                row["contig"],
                str(row["ref_start"]),
                str(row["ref_end"]),
                str(row["intron_index"]),
            ]
            for s in samples:
                fields.append(str(row[s]))
            f.write("\t".join(fields) + "\n")


def write_summary(
    summary_path,
    all_sample_introns,
    ref_introns_total,
    n_introner_removed,
    n_nonintroner,
    n_ref_singletons,
    matrix,
    group1,
    group2,
    folded_afs,
    folded_afs_minor_present,
    folded_afs_minor_absent,
    folded_n,
    unfolded_afs,
    unfolded_afs_loss,
    unfolded_afs_gain,
    unfolded_n,
    n_skipped_outgroup,
    unfolded_1749_afs,
    unfolded_1749_loss,
    unfolded_1749_gain,
    unfolded_1749_n,
    n_skipped_1749,
):
    """Write summary statistics."""
    with open(summary_path, "w") as f:
        f.write("=== Non-introner Intron Polymorphism Summary ===\n\n")

        # Per-sample intron counts
        f.write("Intron counts per sample:\n")
        for sample in sorted(all_sample_introns.keys()):
            n = sum(len(v) for v in all_sample_introns[sample].values())
            n_genes = len(all_sample_introns[sample])
            f.write(f"  {sample}: {n} introns in {n_genes} genes\n")

        f.write(f"\nReference introns (total): {ref_introns_total}\n")
        f.write(f"Introner introns removed: {n_introner_removed}\n")
        f.write(f"Non-introner introns retained: {n_nonintroner}\n")
        f.write(f"Reference-only singletons removed: {n_ref_singletons}\n")

        # Matrix stats
        n_total = len(matrix)
        g1_genotypes_complete = [
            row for row in matrix
            if MISSING not in [row[s] for s in group1]
        ]
        n_complete = len(g1_genotypes_complete)

        n_fixed_present = sum(
            1 for row in g1_genotypes_complete
            if all(row[s] == PRESENT for s in group1)
        )
        n_fixed_absent = sum(
            1 for row in g1_genotypes_complete
            if all(row[s] == ABSENT for s in group1)
        )
        n_polymorphic = n_complete - n_fixed_present - n_fixed_absent

        f.write(f"\nGenotype matrix:\n")
        f.write(f"  Total non-introner intron loci: {n_total}\n")
        f.write(f"  Complete data in Group 1 (no missing): {n_complete}\n")
        f.write(f"  Fixed present (all 11 samples): {n_fixed_present}\n")
        f.write(f"  Fixed absent (all 11 samples): {n_fixed_absent}\n")
        f.write(f"  Polymorphic in Group 1: {n_polymorphic}\n")

        # Folded AFS
        f.write(f"\nFolded AFS (Group 1, n={folded_n}):\n")
        max_maf = folded_n // 2
        f.write(f"  {'MAF':<6} {'Total':>6} {'Minor=Pres':>11} {'Minor=Abs':>10}\n")
        for i in range(1, max_maf + 1):
            total = folded_afs.get(i, 0)
            mp = folded_afs_minor_present.get(i, 0)
            ma = folded_afs_minor_absent.get(i, 0)
            f.write(f"  {i:<6} {total:>6} {mp:>11} {ma:>10}\n")
        f.write(f"  {'Total':<6} {sum(folded_afs.values()):>6} "
                f"{sum(folded_afs_minor_present.values()):>11} "
                f"{sum(folded_afs_minor_absent.values()):>10}\n")

        # Unfolded AFS (all Group 2)
        f.write(f"\nUnfolded AFS — Outgroup: Group 2 (n={unfolded_n}):\n")
        f.write(f"  Skipped (outgroup missing/polymorphic): {n_skipped_outgroup}\n")
        f.write(f"  {'DAC':<6} {'Total':>6} {'Loss':>6} {'Gain':>6}\n")
        for i in range(1, unfolded_n):
            total = unfolded_afs.get(i, 0)
            loss = unfolded_afs_loss.get(i, 0)
            gain = unfolded_afs_gain.get(i, 0)
            f.write(f"  {i:<6} {total:>6} {loss:>6} {gain:>6}\n")
        f.write(f"  {'Total':<6} {sum(unfolded_afs.values()):>6} "
                f"{sum(unfolded_afs_loss.values()):>6} {sum(unfolded_afs_gain.values()):>6}\n")

        # Unfolded AFS (RCC1749 only)
        f.write(f"\nUnfolded AFS — Outgroup: RCC1749 only (n={unfolded_1749_n}):\n")
        f.write(f"  Skipped (outgroup missing): {n_skipped_1749}\n")
        f.write(f"  {'DAC':<6} {'Total':>6} {'Loss':>6} {'Gain':>6}\n")
        for i in range(1, unfolded_1749_n):
            total = unfolded_1749_afs.get(i, 0)
            loss = unfolded_1749_loss.get(i, 0)
            gain = unfolded_1749_gain.get(i, 0)
            f.write(f"  {i:<6} {total:>6} {loss:>6} {gain:>6}\n")
        f.write(f"  {'Total':<6} {sum(unfolded_1749_afs.values()):>6} "
                f"{sum(unfolded_1749_loss.values()):>6} {sum(unfolded_1749_gain.values()):>6}\n")


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Non-introner intron presence/absence polymorphism and AFS"
    )
    parser.add_argument("--gtf_dir", required=True, help="Directory with {sample}.gtf files")
    parser.add_argument("--assemblies_dir", required=True, help="Directory with {sample}.vg_paths.fa assemblies")
    parser.add_argument("--introner_loci", required=True, help="BED file of introner loci (CCMP1545 coords, with 100bp flanking)")
    parser.add_argument("--group1", required=True, help="Comma-separated Group 1 sample names")
    parser.add_argument("--group2", required=True, help="Comma-separated Group 2 sample names")
    parser.add_argument("--reference", required=True, help="Reference sample name (must be in Group 1)")
    parser.add_argument("--output_matrix", required=True, help="Output genotype matrix TSV")
    parser.add_argument("--output_afs", required=True, help="Output AFS plot PDF (also writes sibling PNG)")
    parser.add_argument("--output_summary", required=True, help="Output summary statistics")
    parser.add_argument("--coverage_dir", default=None,
                        help="Directory with {sample}.intron_coverage_calls.tsv files for coverage validation")
    parser.add_argument("--gtf_overrides", default=None,
                        help="Comma-separated sample=path overrides for GTF lookup, "
                             "e.g. RCC1749=/path/to/RCC1749.augmented.renamed.gtf")
    parser.add_argument("--include_mating_region", action="store_true",
                        help="Keep non-introner introns in the CCMP1545 scaffold 2 mating-type region")
    parser.add_argument("--keep_reference_singletons", action="store_true",
                        help="Keep reference-only singleton loci instead of filtering them as likely artifacts")
    args = parser.parse_args()

    group1 = args.group1.split(",")
    group2 = args.group2.split(",")
    reference = args.reference
    all_samples = group1 + group2

    if reference not in group1:
        sys.exit(f"Error: reference {reference} must be in group1")

    gtf_dir = Path(args.gtf_dir)
    assemblies_dir = args.assemblies_dir

    # Parse GTF overrides (e.g. RCC1749=/path/to/augmented.gtf)
    gtf_override_map = {}
    if args.gtf_overrides:
        for entry in args.gtf_overrides.split(","):
            sample_name, gtf_path = entry.split("=", 1)
            gtf_override_map[sample_name] = Path(gtf_path)

    # ── Stage 1: Extract introns from all GTFs ──
    print("Stage 1: Extracting introns from GTFs...")
    all_sample_introns = {}
    for sample in all_samples:
        if sample in gtf_override_map:
            gtf_path = gtf_override_map[sample]
        else:
            gtf_path = gtf_dir / f"{sample}.gtf"
        if not gtf_path.exists():
            sys.exit(f"Error: GTF not found: {gtf_path}")
        introns = extract_introns_from_gtf(gtf_path)
        all_sample_introns[sample] = introns
        n_introns = sum(len(v) for v in introns.values())
        print(f"  {sample}: {len(introns)} genes with introns, {n_introns} total introns")

    # ── Stage 2: Extract intron sequences from assemblies ──
    print("\nStage 2: Extracting intron sequences from assemblies...")
    all_intron_seqs = {}
    for sample in all_samples:
        seqs = extract_intron_sequences(all_sample_introns[sample], assemblies_dir, sample)
        all_intron_seqs[sample] = seqs
        n_seqs = sum(len(v) for v in seqs.values())
        print(f"  {sample}: {n_seqs} intron sequences extracted")

    # ── Stage 3: Filter introner introns from reference ──
    print("\nStage 3: Filtering introner introns from reference...")
    ref_introns = all_sample_introns[reference]
    ref_introns_total = sum(len(v) for v in ref_introns.values())
    print(f"  Reference introns: {ref_introns_total}")

    introner_loci = load_introner_loci(args.introner_loci)
    print(f"  Loaded {len(introner_loci)} introner loci (body coords after stripping flanking)")
    nonintroner_introns, n_removed, n_kept = filter_introner_introns(ref_introns, introner_loci)
    print(f"  Removed {n_removed} introner-overlapping introns")
    print(f"  Retained {n_kept} non-introner introns")

    # Filter out genes in the mating-type region unless explicitly retained.
    if args.include_mating_region:
        print("  Keeping mating-type region introns")
    else:
        MATING_CHROM = "CCMP1545#0#scaffold_2"
        MATING_START = 49808
        MATING_END = 1730591
        n_before = sum(len(v) for v in nonintroner_introns.values())
        nonintroner_introns = {
            gene_id: intron_list
            for gene_id, intron_list in nonintroner_introns.items()
            if not any(
                contig == MATING_CHROM and istart <= MATING_END and iend >= MATING_START
                for contig, istart, iend, _ in intron_list
            )
        }
        n_after = sum(len(v) for v in nonintroner_introns.values())
        print(f"  Removed {n_before - n_after} introns in mating-type region")

    # ── Stage 3b: Load coverage calls (if provided) ──
    coverage_calls = None
    if args.coverage_dir:
        print("\nStage 3b: Loading coverage calls...")
        non_ref_samples = [s for s in all_samples if s != reference]
        coverage_calls = load_coverage_calls(args.coverage_dir, non_ref_samples)
        print(f"  Total coverage calls loaded: {len(coverage_calls)}")

    # ── Stage 4: minimap2 ortholog matching ──
    print("\nStage 4: Running minimap2 ortholog matching...")
    ref_sequences = all_intron_seqs[reference]
    genotype_results = run_minimap2_matching(
        nonintroner_introns, ref_sequences, all_intron_seqs,
        all_sample_introns, group1, group2, reference,
        coverage_calls=coverage_calls,
    )

    # ── Stage 5: Build matrix, compute AFS, plot ──
    print("\nStage 5: Building matrix and computing AFS...")
    matrix = build_genotype_matrix(nonintroner_introns, genotype_results, group1, group2)
    print(f"  Matrix rows (before filtering): {len(matrix)}")

    # Filter reference singletons: introns where CCMP1545 is the only Group 1 sample
    # called PRESENT (others are all ABSENT or MISSING within Group 1).
    # These are likely annotation artifacts in the reference assembly.
    # We check Group 1 only because the outgroup state shouldn't save a ref artifact.
    non_ref_g1 = [s for s in group1 if s != reference]
    n_ref_singletons = 0
    if args.keep_reference_singletons:
        print("  Keeping reference-only singletons")
    else:
        filtered_matrix = []
        for row in matrix:
            if (row[reference] == PRESENT
                    and not any(row[s] == PRESENT for s in non_ref_g1)):
                n_ref_singletons += 1
            else:
                filtered_matrix.append(row)
        matrix = filtered_matrix
        print(f"  Removed {n_ref_singletons} reference-only singletons (likely annotation artifacts)")
    print(f"  Matrix rows: {len(matrix)}")

    folded_afs, folded_afs_minor_present, folded_afs_minor_absent, folded_n = compute_folded_afs(matrix, group1)
    unfolded_afs, unfolded_afs_loss, unfolded_afs_gain, unfolded_n, n_skipped_outgroup = compute_unfolded_afs(matrix, group1, group2)

    # Also compute unfolded AFS using only RCC1749 as outgroup
    rcc1749_outgroup = ["RCC1749"]
    unfolded_1749_afs, unfolded_1749_loss, unfolded_1749_gain, unfolded_1749_n, n_skipped_1749 = compute_unfolded_afs(matrix, group1, rcc1749_outgroup)

    total_folded_poly = sum(folded_afs.values())
    total_unfolded_poly = sum(unfolded_afs.values())
    total_loss = sum(unfolded_afs_loss.values())
    total_gain = sum(unfolded_afs_gain.values())
    total_1749_poly = sum(unfolded_1749_afs.values())
    total_1749_loss = sum(unfolded_1749_loss.values())
    total_1749_gain = sum(unfolded_1749_gain.values())
    print(f"  Folded AFS: {total_folded_poly} polymorphic loci")
    print(f"  Unfolded AFS (Group 2): {total_unfolded_poly} polarizable ({total_loss} loss, {total_gain} gain)")
    print(f"  Unfolded AFS (RCC1749): {total_1749_poly} polarizable ({total_1749_loss} loss, {total_1749_gain} gain)")
    print(f"  Skipped (Group 2 outgroup): {n_skipped_outgroup}")
    print(f"  Skipped (RCC1749 outgroup): {n_skipped_1749}")

    # Write outputs
    Path(args.output_matrix).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output_afs).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output_summary).parent.mkdir(parents=True, exist_ok=True)

    write_matrix(matrix, group1, group2, args.output_matrix)
    print(f"\n  Matrix written to: {args.output_matrix}")

    plot_afs(folded_afs_minor_present, folded_afs_minor_absent, folded_n,
             unfolded_afs_loss, unfolded_afs_gain, unfolded_n,
             unfolded_1749_loss, unfolded_1749_gain, unfolded_1749_n,
             args.output_afs)
    print(f"  AFS plot written to: {args.output_afs}")

    write_summary(
        args.output_summary,
        all_sample_introns,
        ref_introns_total,
        n_removed,
        n_kept,
        n_ref_singletons,
        matrix,
        group1,
        group2,
        folded_afs,
        folded_afs_minor_present,
        folded_afs_minor_absent,
        folded_n,
        unfolded_afs,
        unfolded_afs_loss,
        unfolded_afs_gain,
        unfolded_n,
        n_skipped_outgroup,
        unfolded_1749_afs,
        unfolded_1749_loss,
        unfolded_1749_gain,
        unfolded_1749_n,
        n_skipped_1749,
    )
    print(f"  Summary written to: {args.output_summary}")
    print("\nDone.")


if __name__ == "__main__":
    main()
