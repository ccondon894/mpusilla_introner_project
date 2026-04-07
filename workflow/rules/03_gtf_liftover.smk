# ============================================================
# 03_gtf_liftover.smk - GTF Annotation Liftover
# ============================================================
#
# Lifts over CCMP1545 reference annotations to all sample genomes
# using miniprot protein-to-genome alignment.
#
# Steps:
# 1. Extract reference proteins from CCMP1545 GTF
# 2. Align proteins to each sample genome with miniprot
# 3. Resolve paralog mappings
# 4. Filter low-quality mappings
# 5. Rename gene IDs to match reference
#
# NOTE: The annotate_introners step (update_gtf_with_introners.py)
# has been REMOVED to eliminate circular dependency with genotype matrix.
#
# Adapted from: /scratch1/chris/introner-expression-analysis/rules/make_sample_gtfs.smk
#
# ============================================================

import os
from pathlib import Path

# ============================================================
# CONFIGURATION
# ============================================================

# Reference files
REF_GTF = Path(config["paths"]["references"].get("ccmp1545_gtf",
    "/scratch1/chris/introner-expression-analysis/gene_annotations/MpusillaCCMP1545_228_v3.0.gene_exons.gtf"))
REF_GENOME = Path(config["paths"]["references"]["ccmp1545"])

# Directories
PROTEIN_DIR = RESULTS / "proteins"
STATS_DIR = RESULTS / "annotations" / "stats"

# ============================================================
# RULES
# ============================================================

rule extract_reference_proteins:
    """
    Extract protein sequences from CCMP1545 reference annotation
    using AGAT toolkit
    """
    input:
        gtf = REF_GTF,
        genome = REF_GENOME
    output:
        proteins = PROTEIN_DIR / "reference_proteins.fa"
    shell:
        """
        mkdir -p $(dirname {output.proteins})
        agat_sp_extract_sequences.pl \
            --gff {input.gtf} \
            --fasta {input.genome} \
            --type CDS \
            --merge \
            --protein \
            --output {output.proteins}
        """


rule miniprot_align:
    """
    Align reference proteins to target genome with miniprot
    """
    input:
        proteins = PROTEIN_DIR / "reference_proteins.fa",
        target_genome = ASSEMBLIES_DIR / "{sample}.vg_paths.fa"
    output:
        gtf = ANNOTATIONS_DIR / "{sample}.miniprot.gtf"
    threads: 8
    shell:
        """
        miniprot -t {threads} --gtf --aln {input.target_genome} {input.proteins} > {output.gtf}
        """


rule resolve_paralog_mappings:
    """
    Resolve cases where proteins map to multiple locations
    Keeps the best mapping based on alignment scores
    """
    input:
        miniprot_gtf = ANNOTATIONS_DIR / "{sample}.miniprot.gtf",
        ref_gtf = REF_GTF,
        target_genome = ASSEMBLIES_DIR / "{sample}.vg_paths.fa",
        ref_genome = REF_GENOME
    output:
        resolved_gtf = ANNOTATIONS_DIR / "{sample}.miniprot.resolved.gtf"
    shell:
        """
        python {PROJECT_ROOT}/scripts/genotyping/resolve_paralog_mappings.py \
            --ref-gtf {input.ref_gtf} \
            --ref-fasta {input.ref_genome} \
            --miniprot-gtf {input.miniprot_gtf} \
            --target-fasta {input.target_genome} \
            --output-gtf {output.resolved_gtf}
        """


rule filter_gtf:
    """
    Filter low-quality mappings from resolved GTF
    """
    input:
        gtf = ANNOTATIONS_DIR / "{sample}.miniprot.resolved.gtf",
        ref_gtf = REF_GTF
    output:
        filtered_gtf = ANNOTATIONS_DIR / "{sample}.filtered.gtf"
    shell:
        """
        python {PROJECT_ROOT}/scripts/genotyping/filter_protein_gtf.py \
            --input {input.gtf} \
            --output {output.filtered_gtf}
        """


rule rename_gene_ids:
    """
    Rename gene IDs to match reference annotation format
    """
    input:
        gtf = ANNOTATIONS_DIR / "{sample}.filtered.gtf",
        ref_gtf = REF_GTF
    output:
        gtf = ANNOTATIONS_DIR / "{sample}.gtf"
    shell:
        """
        python {PROJECT_ROOT}/scripts/genotyping/rename_gene_ids.py \
            --reference {input.ref_gtf} \
            --query {input.gtf} \
            --output {output.gtf}
        """


rule compute_gtf_stats:
    """
    Compute statistics comparing query GTF to reference
    """
    input:
        query_gtf = ANNOTATIONS_DIR / "{sample}.gtf",
        ref_gtf = REF_GTF
    output:
        basic = STATS_DIR / "{sample}_basic_stats.tsv",
        ref_chrom = STATS_DIR / "{sample}_ref_chrom_dist.tsv",
        query_chrom = STATS_DIR / "{sample}_query_chrom_dist.tsv",
        overlap = STATS_DIR / "{sample}_overlap_stats.tsv"
    params:
        output_prefix = lambda wildcards: str(STATS_DIR / wildcards.sample)
    shell:
        """
        mkdir -p {STATS_DIR}
        python {PROJECT_ROOT}/scripts/genotyping/gtf_stats.py \
            --reference {input.ref_gtf} \
            --query {input.query_gtf} \
            --output-prefix {params.output_prefix}
        """


rule all_gtf_liftover:
    """
    Target: Generate all GTF annotations
    """
    input:
        expand(ANNOTATIONS_DIR / "{sample}.gtf", sample=ALL_SAMPLES),
        expand(STATS_DIR / "{sample}_basic_stats.tsv", sample=ALL_SAMPLES)
