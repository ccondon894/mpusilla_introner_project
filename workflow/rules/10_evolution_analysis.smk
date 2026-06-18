# ============================================================
# 10_evolution_analysis.smk - Evolution & Diversity Analysis
# ============================================================
#
# Analyzes evolutionary patterns around introner insertion sites
# by examining flanking sequence diversity.
#
# Contains two analysis modes:
# 1. All-Samples Analysis: Compares fixation patterns between Group1/Group2
# 2. Group1 Analysis: Frequency-based polymorphism analysis within Group1
#
# Steps:
# - Build consensus sequences from read alignments
# - Classify orthologs by presence/absence patterns
# - Align flanking regions with MAFFT
# - Calculate nucleotide diversity (π) and between-group divergence (dxy)
# - Generate visualization plots
#
# Adapted from:
# - /scratch1/chris/introner-genotyping-pipeline/rules/all_samples_evolution_analysis.smk
# - /scratch1/chris/introner-genotyping-pipeline/rules/group1_evolution_analysis.smk
#
# ============================================================

import os
import glob
import json
from pathlib import Path

# ============================================================
# CONFIGURATION
# ============================================================

# Output directories
EVOLUTION_DIR = RESULTS / "evolution"
CONSENSUS_DIR = EVOLUTION_DIR / "consensus_sequences"
ALIGNMENT_DIR = EVOLUTION_DIR / "alignments"
DIVERSITY_DIR = EVOLUTION_DIR / "diversity_metrics"
EVOLUTION_PLOTS_DIR = EVOLUTION_DIR / "plots"
EVOLUTION_LOG_DIR = EVOLUTION_DIR / "logs"
INTRONER_BODY_DECAY_DIR = EVOLUTION_DIR / "introner_body_decay"
INTRONER_BODY_CONSENSUS_DIR = EVOLUTION_DIR / "introner_body_consensus"

# Local aliases: later-loaded rule files (30_expression.smk,
# 31_isoform_analysis.smk) redefine ALIGNMENT_DIR and DIVERSITY_DIR at
# global scope. We capture the evolution paths here so lambdas and shell
# directives within this file reference the correct directory regardless
# of load order.
_EVO_ALIGNMENT_DIR = ALIGNMENT_DIR
_EVO_DIVERSITY_DIR = DIVERSITY_DIR
INTRONER_BODY_DECAY_ALIGNMENT_DIR = _EVO_ALIGNMENT_DIR / "introner_body_decay"
INTRONER_BODY_CONSENSUS_QC = INTRONER_BODY_CONSENSUS_DIR / "introner_body_consensus.qc.tsv"
EVO_NON_REF_SAMPLES = [sample for sample in ALL_SAMPLES if sample != REFERENCE]

# Input directories from previous steps
COVERAGE_BAM_DIR = GENOTYPING_DIR / "coverage" / "bams"

# Flanking sequence lengths to analyze
FLANK_LENGTHS = config["params"]["flanks"].get("analysis_lengths", [100, 200])
SIDES = ["left", "right"]

# Fixation categories for all-samples analysis
FIXATION_CATEGORIES = ["group1_fixed_group2_absent", "group1_absent_group2_fixed", "group1_fixed_group2_fixed"]

# Frequency categories for Group1 analysis (1-10 out of 11 samples)
FREQUENCY_CATEGORIES = list(range(1, 11))

# Wildcard constraints
wildcard_constraints:
    flank_length = "|".join(map(str, FLANK_LENGTHS)),
    category = "|".join(FIXATION_CATEGORIES),
    freq = "|".join(map(str, FREQUENCY_CATEGORIES)),
    side = "left|right",
    group = "group1|group2",
    state = "present|absent"


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def _collect_all_samples_alignment_files(classification_file, flank_length):
    """Return alignment files whose source .fa was actually written.

    classify_all_samples_orthologs.py only writes a per-group .fa file when
    at least one sequence could be extracted (e.g. absent samples without
    valid coords contribute nothing, so the file may not be created). We
    filter to existing .fa sources here so Snakemake doesn't request
    non-buildable .mafft.fa files.
    """
    with open(classification_file, "r") as f:
        classification = json.load(f)

    alignment_dir_bp = _EVO_ALIGNMENT_DIR / f"{flank_length}bp"
    alignment_files = []

    for oid, info in classification.items():
        category = info["category"]
        for side in SIDES:
            g1_src = alignment_dir_bp / f"{oid}.{category}.group1.{side}_flank_{flank_length}bp.fa"
            g2_src = alignment_dir_bp / f"{oid}.{category}.group2.{side}_flank_{flank_length}bp.fa"

            if g1_src.exists():
                alignment_files.append(
                    alignment_dir_bp / f"{oid}.{category}.group1.{side}_flank_{flank_length}bp.mafft.fa"
                )
            if g2_src.exists():
                alignment_files.append(
                    alignment_dir_bp / f"{oid}.{category}.group2.{side}_flank_{flank_length}bp.mafft.fa"
                )
            if g1_src.exists() and g2_src.exists():
                alignment_files.append(
                    alignment_dir_bp / f"{oid}.{category}.combined.{side}_flank_{flank_length}bp.mafft.fa"
                )

    return alignment_files


def get_all_samples_alignment_files(wildcards):
    """Get all alignment files for all-samples analysis"""
    classification_file = _EVO_ALIGNMENT_DIR / f"all_samples_classification_{wildcards.flank_length}bp.json"
    if exists(classification_file):
        return _collect_all_samples_alignment_files(classification_file, wildcards.flank_length)

    # Defer evaluation until checkpoint completes - this prevents file I/O during DAG construction
    checkpoints.classify_all_samples_orthologs.get()
    if not classification_file.exists():
        return []
    return _collect_all_samples_alignment_files(classification_file, wildcards.flank_length)


def _collect_group1_alignment_files(classification_file, flank_length):
    """Return Group1 alignment files whose source .fa was actually written."""
    with open(classification_file, "r") as f:
        classification = json.load(f)

    alignment_dir_bp = _EVO_ALIGNMENT_DIR / f"{flank_length}bp"
    alignment_files = []

    for oid, info in classification.items():
        freq = info["frequency"]
        for side in SIDES:
            present_src = alignment_dir_bp / f"{oid}.freq_{freq}.present.{side}_flank_{flank_length}bp.fa"
            absent_src = alignment_dir_bp / f"{oid}.freq_{freq}.absent.{side}_flank_{flank_length}bp.fa"

            if present_src.exists():
                alignment_files.append(
                    alignment_dir_bp / f"{oid}.freq_{freq}.present.{side}_flank_{flank_length}bp.mafft.fa"
                )
            if absent_src.exists():
                alignment_files.append(
                    alignment_dir_bp / f"{oid}.freq_{freq}.absent.{side}_flank_{flank_length}bp.mafft.fa"
                )
            if present_src.exists() and absent_src.exists():
                alignment_files.append(
                    alignment_dir_bp / f"{oid}.freq_{freq}.combined.{side}_flank_{flank_length}bp.mafft.fa"
                )

    return alignment_files


def get_group1_alignment_files(wildcards):
    """Get all alignment files for Group1 analysis"""
    classification_file = _EVO_ALIGNMENT_DIR / f"group1_classification_{wildcards.flank_length}bp.json"

    if exists(classification_file):
        return _collect_group1_alignment_files(classification_file, wildcards.flank_length)

    # Defer evaluation until checkpoint completes
    checkpoints.classify_group1_orthologs.get()
    if not classification_file.exists():
        return []
    return _collect_group1_alignment_files(classification_file, wildcards.flank_length)


# ============================================================
# CONSENSUS SEQUENCE BUILDING
# ============================================================

rule make_evolution_bed_files:
    """
    Create BED files for flanking regions around introner loci.

    Extracts coordinates for left and right flanking regions at
    specified lengths (100bp, 200bp) for each sample.
    """
    input:
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv"
    output:
        left_bed = CONSENSUS_DIR / "{sample}.loci.left_flank_{flank_length}bp.bed",
        right_bed = CONSENSUS_DIR / "{sample}.loci.right_flank_{flank_length}bp.bed"
    shell:
        """
        mkdir -p {CONSENSUS_DIR}
        python {PROJECT_ROOT}/scripts/evolution/make_all_samples_loci_flank_beds.py \
            {input.matrix} {wildcards.sample} {wildcards.flank_length} \
            {output.left_bed} {output.right_bed}
        """


rule evolution_mpileup:
    """
    Generate mpileup for flanking regions to build consensus.

    Uses bcftools mpileup with quality filters:
    -d 100: max depth
    -q 30: min mapping quality
    -Q 20: min base quality
    """
    input:
        bam = COVERAGE_BAM_DIR / "{sample}.sorted.bam",
        bed = CONSENSUS_DIR / "{sample}.loci.{side}_flank_{flank_length}bp.bed",
        fa = ASSEMBLIES_DIR / "{sample}.vg_paths.fa"
    output:
        mpileup = CONSENSUS_DIR / "{sample}.loci.{side}_flank_{flank_length}bp.mpileup.bcf"
    params:
        temp_dir = lambda wildcards: f"temp_{wildcards.sample}_{wildcards.side}_{wildcards.flank_length}bp_sort"
    shell:
        """
        mkdir -p {params.temp_dir}

        bcftools mpileup \
            -d 100 -q 30 -Q 20 -A \
            -f {input.fa} -R {input.bed} \
            -Ou {input.bam} | \
        bcftools sort -T {params.temp_dir} -Ob -o {output.mpileup}

        rm -rf {params.temp_dir}
        """


rule evolution_call_variants:
    """
    Call variants from mpileup for consensus building.

    Uses haploid mode (--ploidy 1) and filters for depth >= 10.
    """
    input:
        mpileup = CONSENSUS_DIR / "{sample}.loci.{side}_flank_{flank_length}bp.mpileup.bcf",
        fa = ASSEMBLIES_DIR / "{sample}.vg_paths.fa"
    output:
        vcf = CONSENSUS_DIR / "{sample}.loci.{side}_flank_{flank_length}bp.vcf.gz",
        filt_vcf = CONSENSUS_DIR / "{sample}.loci.{side}_flank_{flank_length}bp.filtered.vcf.gz"
    shell:
        """
        bcftools call -m --ploidy 1 -Ou {input.mpileup} | \
        bcftools norm -f {input.fa} -m +both -Oz -o {output.vcf}

        bcftools filter -i 'DP>=10' {output.vcf} -o {output.filt_vcf}
        tabix -p vcf {output.filt_vcf}
        """


rule build_sample_consensus:
    """
    Build consensus sequences for flanking regions.

    Applies variants from filtered VCF to reference to create
    sample-specific consensus sequences for each locus.
    """
    input:
        filt_vcf = CONSENSUS_DIR / "{sample}.loci.{side}_flank_{flank_length}bp.filtered.vcf.gz",
        bed = CONSENSUS_DIR / "{sample}.loci.{side}_flank_{flank_length}bp.bed",
        fa = ASSEMBLIES_DIR / "{sample}.vg_paths.fa"
    output:
        cons = CONSENSUS_DIR / "{sample}.loci.{side}_flank_{flank_length}bp.consensus.fa"
    params:
        cons_genome = lambda wildcards: str(CONSENSUS_DIR / f"{wildcards.sample}.{wildcards.side}_{wildcards.flank_length}bp.consensus.fa"),
        loci = lambda wildcards: str(CONSENSUS_DIR / f"{wildcards.sample}.loci.{wildcards.side}_flank_{wildcards.flank_length}bp.txt")
    shell:
        """
        bcftools consensus -a N -f {input.fa} -o {params.cons_genome} {input.filt_vcf}
        samtools faidx {params.cons_genome}
        awk '{{print $1":"$2"-"$3}}' {input.bed} > {params.loci}
        samtools faidx {params.cons_genome} $(cat {params.loci}) > {output.cons}
        rm -f {params.cons_genome} {params.cons_genome}.fai {params.loci}
        """


rule build_ref_consensus:
    """
    Extract reference consensus sequences (no variant calling needed).
    """
    input:
        bed = CONSENSUS_DIR / f"{REFERENCE}.loci.{{side}}_flank_{{flank_length}}bp.bed",
        fa = ASSEMBLIES_DIR / f"{REFERENCE}.vg_paths.fa"
    output:
        fa = CONSENSUS_DIR / f"{REFERENCE}.loci.{{side}}_flank_{{flank_length}}bp.consensus.fa"
    shell:
        """
        bedtools getfasta -fi {input.fa} -bed {input.bed} -fo {output.fa} -name
        sed -i 's/ortholog_id_[0-9]\\{{4\\}}:://g' {output.fa}
        """


rule make_introner_body_bed:
    """
    Create BED files for present introner bodies, excluding 100 bp matrix flanks.
    """
    input:
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv"
    output:
        bed = INTRONER_BODY_CONSENSUS_DIR / "{sample}.introner_body.bed"
    shell:
        """
        mkdir -p {INTRONER_BODY_CONSENSUS_DIR}
        python {PROJECT_ROOT}/scripts/evolution/make_introner_body_bed.py \
            {input.matrix} {wildcards.sample} {output.bed} \
            --flank-length 100
        """


rule introner_body_mpileup:
    """
    Generate mpileup for non-reference introner bodies.
    """
    input:
        bam = COVERAGE_BAM_DIR / "{sample}.sorted.bam",
        bed = INTRONER_BODY_CONSENSUS_DIR / "{sample}.introner_body.bed",
        fa = ASSEMBLIES_DIR / "{sample}.vg_paths.fa"
    output:
        mpileup = INTRONER_BODY_CONSENSUS_DIR / "{sample}.introner_body.mpileup.bcf"
    params:
        temp_dir = lambda wc: f"/scratch1/chris/tmp/introner_body_{wc.sample}_sort"
    wildcard_constraints:
        sample = "|".join(EVO_NON_REF_SAMPLES)
    shell:
        """
        mkdir -p {params.temp_dir}
        bcftools mpileup \
            -d 100 -q 30 -Q 20 -A \
            -f {input.fa} -R {input.bed} \
            -Ou {input.bam} | \
        bcftools sort -T {params.temp_dir} -Ob -o {output.mpileup}
        rm -rf {params.temp_dir}
        """


rule introner_body_call_variants:
    """
    Call haploid variants for non-reference introner-body consensus.
    """
    input:
        mpileup = INTRONER_BODY_CONSENSUS_DIR / "{sample}.introner_body.mpileup.bcf",
        fa = ASSEMBLIES_DIR / "{sample}.vg_paths.fa"
    output:
        vcf = INTRONER_BODY_CONSENSUS_DIR / "{sample}.introner_body.vcf.gz",
        filt_vcf = INTRONER_BODY_CONSENSUS_DIR / "{sample}.introner_body.filtered.vcf.gz"
    wildcard_constraints:
        sample = "|".join(EVO_NON_REF_SAMPLES)
    shell:
        """
        bcftools call -m --ploidy 1 -Ou {input.mpileup} | \
        bcftools norm -f {input.fa} -m +both -Oz -o {output.vcf}
        bcftools filter -i 'DP>=10' {output.vcf} -o {output.filt_vcf}
        tabix -p vcf {output.filt_vcf}
        """


rule build_introner_body_consensus:
    """
    Build non-reference introner-body consensus FASTA from filtered variants.
    """
    input:
        filt_vcf = INTRONER_BODY_CONSENSUS_DIR / "{sample}.introner_body.filtered.vcf.gz",
        bed = INTRONER_BODY_CONSENSUS_DIR / "{sample}.introner_body.bed",
        fa = ASSEMBLIES_DIR / "{sample}.vg_paths.fa"
    output:
        cons = INTRONER_BODY_CONSENSUS_DIR / "{sample}.introner_body.consensus.fa"
    params:
        cons_genome = lambda wc: str(INTRONER_BODY_CONSENSUS_DIR / f"{wc.sample}.introner_body.tmp_consensus.fa")
    wildcard_constraints:
        sample = "|".join(EVO_NON_REF_SAMPLES)
    shell:
        """
        bcftools consensus -a N -f {input.fa} -o {params.cons_genome} {input.filt_vcf}
        bedtools getfasta -fi {params.cons_genome} -bed {input.bed} -fo {output.cons} -nameOnly
        rm -f {params.cons_genome} {params.cons_genome}.fai
        """


rule build_ref_introner_body_consensus:
    """
    Extract reference introner-body sequences directly from the reference assembly.
    """
    input:
        bed = INTRONER_BODY_CONSENSUS_DIR / f"{REFERENCE}.introner_body.bed",
        fa = ASSEMBLIES_DIR / f"{REFERENCE}.vg_paths.fa"
    output:
        cons = INTRONER_BODY_CONSENSUS_DIR / f"{REFERENCE}.introner_body.consensus.fa"
    shell:
        """
        bedtools getfasta -fi {input.fa} -bed {input.bed} -fo {output.cons} -nameOnly
        """


rule qc_introner_body_consensus:
    """
    Apply strict repetitive-region QC to non-reference introner-body consensus.
    """
    input:
        bed = INTRONER_BODY_CONSENSUS_DIR / "{sample}.introner_body.bed",
        cons = INTRONER_BODY_CONSENSUS_DIR / "{sample}.introner_body.consensus.fa",
        bam = COVERAGE_BAM_DIR / "{sample}.sorted.bam",
        fa = ASSEMBLIES_DIR / "{sample}.vg_paths.fa"
    output:
        qc = INTRONER_BODY_CONSENSUS_DIR / "{sample}.introner_body.qc.tsv"
    wildcard_constraints:
        sample = "|".join(EVO_NON_REF_SAMPLES)
    shell:
        """
        python {PROJECT_ROOT}/scripts/evolution/qc_introner_body_consensus.py \
            --sample {wildcards.sample} \
            --bed {input.bed} \
            --consensus-fasta {input.cons} \
            --bam {input.bam} \
            --reference-fasta {input.fa} \
            --output {output.qc} \
            --reference-sample {REFERENCE}
        """


rule qc_ref_introner_body_consensus:
    """
    Mark reference introner-body consensus sequences as reference-derived QC pass.
    """
    input:
        bed = INTRONER_BODY_CONSENSUS_DIR / f"{REFERENCE}.introner_body.bed",
        cons = INTRONER_BODY_CONSENSUS_DIR / f"{REFERENCE}.introner_body.consensus.fa",
        fa = ASSEMBLIES_DIR / f"{REFERENCE}.vg_paths.fa"
    output:
        qc = INTRONER_BODY_CONSENSUS_DIR / f"{REFERENCE}.introner_body.qc.tsv"
    shell:
        """
        python {PROJECT_ROOT}/scripts/evolution/qc_introner_body_consensus.py \
            --sample {REFERENCE} \
            --bed {input.bed} \
            --consensus-fasta {input.cons} \
            --reference-fasta {input.fa} \
            --output {output.qc} \
            --reference-sample {REFERENCE}
        """


rule aggregate_introner_body_consensus_qc:
    """
    Combine per-sample introner-body QC tables.
    """
    input:
        qcs = expand(INTRONER_BODY_CONSENSUS_DIR / "{sample}.introner_body.qc.tsv", sample=ALL_SAMPLES)
    output:
        qc = INTRONER_BODY_CONSENSUS_QC
    shell:
        """
        awk 'FNR == 1 && NR != 1 {{next}} {{print}}' {input.qcs} > {output.qc}
        """





# ============================================================
# ALL-SAMPLES ANALYSIS (Group1 vs Group2 Fixation Patterns)
# ============================================================

checkpoint classify_all_samples_orthologs:
    """
    Classify orthologs by fixation status between Group1 and Group2.

    Categories:
    - group1_fixed_group2_absent: Fixed in Group1, absent in Group2
    - group1_absent_group2_fixed: Absent in Group1, fixed in Group2
    - group1_fixed_group2_fixed: Fixed in both groups

    Creates FASTA files for each ortholog grouped by category.
    """
    input:
        fastas = expand(
            CONSENSUS_DIR / "{sample}.loci.{side}_flank_{flank_length}bp.consensus.fa",
            sample=ALL_SAMPLES, side=SIDES, flank_length=FLANK_LENGTHS
        ),
        gt_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv"
    output:
        classification_100bp = ALIGNMENT_DIR / "all_samples_classification_100bp.json",
        classification_200bp = ALIGNMENT_DIR / "all_samples_classification_200bp.json"
    params:
        consensus_dir = CONSENSUS_DIR,
        alignment_dir = ALIGNMENT_DIR,
        group1_str = ",".join(GROUP1_SAMPLES),
        group2_str = ",".join(GROUP2_SAMPLES)
    shell:
        """
        mkdir -p {params.alignment_dir}/100bp
        mkdir -p {params.alignment_dir}/200bp

        python {PROJECT_ROOT}/scripts/evolution/classify_all_samples_orthologs.py \
            {input.gt_matrix} {params.consensus_dir} {params.alignment_dir} \
            --group1 {params.group1_str} \
            --group2 {params.group2_str}
        """


rule align_all_samples_groups:
    """
    Align flanking sequences within each group using MAFFT.
    """
    input:
        fasta = ALIGNMENT_DIR / "{flank_length}bp" / "{ortholog_id}.{category}.{group}.{side}_flank_{flank_length}bp.fa"
    output:
        aligned = ALIGNMENT_DIR / "{flank_length}bp" / "{ortholog_id}.{category}.{group}.{side}_flank_{flank_length}bp.mafft.fa"
    wildcard_constraints:
        group = "group1|group2"
    shell:
        """
        mafft --adjustdirection --maxiterate 1000 --globalpair --quiet {input.fasta} > {output.aligned}
        """


rule align_all_samples_combined:
    """
    Create combined alignments for between-group divergence (dxy) calculation.
    """
    input:
        group1 = ALIGNMENT_DIR / "{flank_length}bp" / "{ortholog_id}.{category}.group1.{side}_flank_{flank_length}bp.fa",
        group2 = ALIGNMENT_DIR / "{flank_length}bp" / "{ortholog_id}.{category}.group2.{side}_flank_{flank_length}bp.fa"
    output:
        combined = ALIGNMENT_DIR / "{flank_length}bp" / "{ortholog_id}.{category}.combined.{side}_flank_{flank_length}bp.fa",
        aligned = ALIGNMENT_DIR / "{flank_length}bp" / "{ortholog_id}.{category}.combined.{side}_flank_{flank_length}bp.mafft.fa"
    shell:
        """
        cat {input.group1} {input.group2} > {output.combined}
        mafft --adjustdirection --maxiterate 1000 --globalpair --quiet {output.combined} > {output.aligned}
        """


rule calculate_all_samples_diversity:
    """
    Calculate diversity metrics (π, dxy) for all-samples analysis.

    Computes:
    - Within-group nucleotide diversity (π) for Group1 and Group2
    - Between-group divergence (dxy)
    """
    input:
        alignment_files = get_all_samples_alignment_files
    output:
        metrics = DIVERSITY_DIR / "all_samples_diversity_metrics_{flank_length}bp.tsv"
    params:
        alignment_dir = _EVO_ALIGNMENT_DIR,
        diversity_dir = _EVO_DIVERSITY_DIR,
        classification = lambda wildcards: str(_EVO_ALIGNMENT_DIR / f"all_samples_classification_{wildcards.flank_length}bp.json")
    shell:
        """
        mkdir -p {params.diversity_dir}
        python {PROJECT_ROOT}/scripts/evolution/calculate_all_samples_diversity_metrics_ros.py \
            --alignment_dir {params.alignment_dir} \
            --classification {params.classification} \
            --flank_length {wildcards.flank_length} \
            --output {output.metrics}
        """


rule calculate_shared_introner_body_dxy:
    """
    Calculate introner-body Dxy for fixed shared loci using active all-samples
    classification classes.
    """
    input:
        classification = _EVO_ALIGNMENT_DIR / "all_samples_classification_{flank_length}bp.json",
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        consensus = expand(INTRONER_BODY_CONSENSUS_DIR / "{sample}.introner_body.consensus.fa", sample=ALL_SAMPLES),
        qc = INTRONER_BODY_CONSENSUS_QC
    output:
        body_dxy = _EVO_DIVERSITY_DIR / "shared_introner_body_dxy_{flank_length}bp.tsv"
    params:
        consensus_dir = INTRONER_BODY_CONSENSUS_DIR,
        diversity_dir = _EVO_DIVERSITY_DIR,
        alignment_dir = _EVO_ALIGNMENT_DIR / "introner_body",
        group1 = ' '.join(GROUP1_SAMPLES),
        group2 = ' '.join(GROUP2_SAMPLES)
    shell:
        """
        mkdir -p {params.alignment_dir} {params.diversity_dir}
        python {PROJECT_ROOT}/scripts/evolution/calculate_shared_introner_body_dxy.py \
            --genotype-matrix {input.genotype_matrix} \
            --classification {input.classification} \
            --consensus-dir {params.consensus_dir} \
            --qc-table {input.qc} \
            --alignment-dir {params.alignment_dir} \
            --output {output.body_dxy} \
            --group1-samples {params.group1} \
            --group2-samples {params.group2}
        """


rule plot_all_samples_analysis:
    """
    Generate the combined all-samples diversity plot.
    """
    input:
        metrics = DIVERSITY_DIR / "all_samples_diversity_metrics_{flank_length}bp.tsv",
        body_dxy = _EVO_DIVERSITY_DIR / "shared_introner_body_dxy_{flank_length}bp.tsv",
        body_decay = INTRONER_BODY_DECAY_DIR / "introner_body_decay.per_locus.tsv"
    output:
        box_plot = EVOLUTION_PLOTS_DIR / "all_samples_box_plots_{flank_length}bp.png"
    shell:
        """
        mkdir -p {EVOLUTION_PLOTS_DIR}
        python {PROJECT_ROOT}/scripts/evolution/plot_all_samples_analysis_boxplot.py \
            --input {input.metrics} \
            --body-dxy {input.body_dxy} \
            --body-decay {input.body_decay} \
            --output {output.box_plot} \
            --flank_length {wildcards.flank_length}
        """


# ============================================================
# GROUP1 ANALYSIS (Frequency-Based Polymorphism)
# ============================================================

checkpoint classify_group1_orthologs:
    """
    Classify Group1 orthologs by frequency (1-10 out of 11 samples).

    For polymorphic introners within Group1, creates FASTA files
    separating samples where introner is present vs absent.
    """
    input:
        fastas = expand(
            CONSENSUS_DIR / "{sample}.loci.{side}_flank_{flank_length}bp.consensus.fa",
            sample=GROUP1_SAMPLES, side=SIDES, flank_length=FLANK_LENGTHS
        ),
        gt_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv"
    output:
        classification_100bp = ALIGNMENT_DIR / "group1_classification_100bp.json",
        classification_200bp = ALIGNMENT_DIR / "group1_classification_200bp.json"
    params:
        consensus_dir = CONSENSUS_DIR,
        alignment_dir = ALIGNMENT_DIR,
        group1_str = ",".join(GROUP1_SAMPLES)
    shell:
        """
        mkdir -p {params.alignment_dir}/100bp
        mkdir -p {params.alignment_dir}/200bp

        python {PROJECT_ROOT}/scripts/evolution/classify_group1_orthologs.py \
            {input.gt_matrix} {params.consensus_dir} {params.alignment_dir} \
            --group1 {params.group1_str}
        """


rule align_group1_frequency_groups:
    """
    Align flanking sequences for present/absent sample sets.
    """
    input:
        fasta = ALIGNMENT_DIR / "{flank_length}bp" / "{ortholog_id}.freq_{freq}.{state}.{side}_flank_{flank_length}bp.fa"
    output:
        aligned = ALIGNMENT_DIR / "{flank_length}bp" / "{ortholog_id}.freq_{freq}.{state}.{side}_flank_{flank_length}bp.mafft.fa"
    wildcard_constraints:
        state = "present|absent"
    shell:
        """
        mafft --adjustdirection --maxiterate 1000 --globalpair --quiet {input.fasta} > {output.aligned}
        """


rule align_group1_combined:
    """
    Create combined alignments for present vs absent comparison.
    """
    input:
        present = ALIGNMENT_DIR / "{flank_length}bp" / "{ortholog_id}.freq_{freq}.present.{side}_flank_{flank_length}bp.fa",
        absent = ALIGNMENT_DIR / "{flank_length}bp" / "{ortholog_id}.freq_{freq}.absent.{side}_flank_{flank_length}bp.fa"
    output:
        combined = ALIGNMENT_DIR / "{flank_length}bp" / "{ortholog_id}.freq_{freq}.combined.{side}_flank_{flank_length}bp.fa",
        aligned = ALIGNMENT_DIR / "{flank_length}bp" / "{ortholog_id}.freq_{freq}.combined.{side}_flank_{flank_length}bp.mafft.fa"
    shell:
        """
        cat {input.present} {input.absent} > {output.combined}
        mafft --adjustdirection --maxiterate 1000 --globalpair --quiet {output.combined} > {output.aligned}
        """


rule calculate_group1_diversity:
    """
    Calculate diversity metrics for Group1 frequency analysis.

    Computes π separately for present vs absent sample sets
    at each frequency category.
    """
    input:
        alignment_files = get_group1_alignment_files
    output:
        metrics = DIVERSITY_DIR / "group1_diversity_metrics_{flank_length}bp.tsv"
    params:
        alignment_dir = _EVO_ALIGNMENT_DIR,
        diversity_dir = _EVO_DIVERSITY_DIR,
        classification = lambda wildcards: str(_EVO_ALIGNMENT_DIR / f"group1_classification_{wildcards.flank_length}bp.json")
    shell:
        """
        mkdir -p {params.diversity_dir}
        python {PROJECT_ROOT}/scripts/evolution/calculate_group1_diversity_metrics_ros.py \
            --alignment_dir {params.alignment_dir} \
            --classification {params.classification} \
            --flank_length {wildcards.flank_length} \
            --output {output.metrics}
        """


rule plot_group1_frequency_lines:
    """
    Generate frequency-based line plots for Group1 analysis.
    """
    input:
        metrics = DIVERSITY_DIR / "group1_diversity_metrics_{flank_length}bp.tsv"
    output:
        present_plot = EVOLUTION_PLOTS_DIR / "frequency_line_plots_{flank_length}bp_present_spectrum.png",
        absent_plot = EVOLUTION_PLOTS_DIR / "frequency_line_plots_{flank_length}bp_absent_spectrum.png"
    shell:
        """
        mkdir -p {EVOLUTION_PLOTS_DIR}
        python {PROJECT_ROOT}/scripts/evolution/plot_group1_frequency_analysis.py \
            --input {input.metrics} \
            --output {EVOLUTION_PLOTS_DIR} \
            --plot_type frequency_lines \
            --flank_length {wildcards.flank_length}
        """


rule plot_group1_box_plots:
    """
    Generate box plots for Group1 frequency analysis.
    """
    input:
        metrics = DIVERSITY_DIR / "group1_diversity_metrics_{flank_length}bp.tsv"
    output:
        present_box = EVOLUTION_PLOTS_DIR / "frequency_box_plots_{flank_length}bp_present_spectrum.png",
        absent_box = EVOLUTION_PLOTS_DIR / "frequency_box_plots_{flank_length}bp_absent_spectrum.png"
    shell:
        """
        mkdir -p {EVOLUTION_PLOTS_DIR}
        python {PROJECT_ROOT}/scripts/evolution/plot_group1_frequency_analysis.py \
            --input {input.metrics} \
            --output {EVOLUTION_PLOTS_DIR} \
            --plot_type box_plot \
            --flank_length {wildcards.flank_length}
        """


# ============================================================
# INTRONER BODY SEQUENCE DECAY
# ============================================================

rule analyze_introner_body_decay:
    """
    Compare introner-body sequence decay between fixed-present and polymorphic introners.

    The primary analysis is Group 1 fixed-present versus complete-call
    polymorphic introners in the three largest present-introner families.
    Group 2 is emitted as a low-power secondary comparison.
    """
    input:
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        consensus = expand(INTRONER_BODY_CONSENSUS_DIR / "{sample}.introner_body.consensus.fa", sample=ALL_SAMPLES),
        qc = INTRONER_BODY_CONSENSUS_QC
    output:
        per_locus = INTRONER_BODY_DECAY_DIR / "introner_body_decay.per_locus.tsv",
        summary = INTRONER_BODY_DECAY_DIR / "introner_body_decay.summary.tsv",
        tests = INTRONER_BODY_DECAY_DIR / "introner_body_decay.tests.tsv",
        plot = EVOLUTION_PLOTS_DIR / "introner_body_decay.top3_families.png"
    params:
        consensus_dir = INTRONER_BODY_CONSENSUS_DIR,
        alignment_dir = INTRONER_BODY_DECAY_ALIGNMENT_DIR,
        group1 = " ".join(GROUP1_SAMPLES),
        group2 = " ".join(GROUP2_SAMPLES),
        top_families = 3
    log:
        EVOLUTION_LOG_DIR / "introner_body_decay.log"
    shell:
        """
        mkdir -p {INTRONER_BODY_DECAY_DIR} {params.alignment_dir} {EVOLUTION_PLOTS_DIR} $(dirname {log})
        python {PROJECT_ROOT}/scripts/evolution/analyze_introner_body_decay.py \
            --genotype-matrix {input.genotype_matrix} \
            --consensus-dir {params.consensus_dir} \
            --qc-table {input.qc} \
            --alignment-dir {params.alignment_dir} \
            --per-locus {output.per_locus} \
            --summary {output.summary} \
            --tests {output.tests} \
            --plot {output.plot} \
            --group1-samples {params.group1} \
            --group2-samples {params.group2} \
            --top-families {params.top_families} \
            > {log} 2>&1
        """


rule summarize_introner_body_decay_tail:
    """
    Summarize high-pi polymorphic tail loci and frequency-level body decay.
    """
    input:
        per_locus = INTRONER_BODY_DECAY_DIR / "introner_body_decay.per_locus.tsv"
    output:
        tail_loci = INTRONER_BODY_DECAY_DIR / "introner_body_decay.polymorphic_pi_tail.tsv",
        frequency_summary = INTRONER_BODY_DECAY_DIR / "introner_body_decay.frequency_summary.tsv",
        tail_summary = INTRONER_BODY_DECAY_DIR / "introner_body_decay.tail_summary.tsv",
        plot = EVOLUTION_PLOTS_DIR / "introner_body_decay.frequency_tail_diagnostics.png"
    log:
        EVOLUTION_LOG_DIR / "introner_body_decay_tail.log"
    shell:
        """
        mkdir -p {INTRONER_BODY_DECAY_DIR} {EVOLUTION_PLOTS_DIR} $(dirname {log})
        python {PROJECT_ROOT}/scripts/evolution/summarize_introner_body_decay_tail.py \
            --per-locus {input.per_locus} \
            --tail-loci {output.tail_loci} \
            --frequency-summary {output.frequency_summary} \
            --tail-summary {output.tail_summary} \
            --plot {output.plot} \
            > {log} 2>&1
        """


# ============================================================
# TARGET RULES
# ============================================================

rule introner_body_decay:
    """
    Target: introner-body sequence decay analysis for fixed vs polymorphic introners.
    """
    input:
        INTRONER_BODY_DECAY_DIR / "introner_body_decay.per_locus.tsv",
        INTRONER_BODY_DECAY_DIR / "introner_body_decay.summary.tsv",
        INTRONER_BODY_DECAY_DIR / "introner_body_decay.tests.tsv",
        EVOLUTION_PLOTS_DIR / "introner_body_decay.top3_families.png",
        INTRONER_BODY_DECAY_DIR / "introner_body_decay.polymorphic_pi_tail.tsv",
        INTRONER_BODY_DECAY_DIR / "introner_body_decay.frequency_summary.tsv",
        INTRONER_BODY_DECAY_DIR / "introner_body_decay.tail_summary.tsv",
        EVOLUTION_PLOTS_DIR / "introner_body_decay.frequency_tail_diagnostics.png"


rule introner_body_consensus_only:
    """
    Target: strict introner-body consensus sequences and QC table.
    """
    input:
        expand(INTRONER_BODY_CONSENSUS_DIR / "{sample}.introner_body.consensus.fa", sample=ALL_SAMPLES),
        INTRONER_BODY_CONSENSUS_QC


rule all_evolution_analysis:
    """
    Target: Complete evolution analysis (both all-samples and Group1)
    """
    input:
        # All-samples analysis outputs
        expand(DIVERSITY_DIR / "all_samples_diversity_metrics_{flank_length}bp.tsv",
               flank_length=FLANK_LENGTHS),
        expand(EVOLUTION_PLOTS_DIR / "all_samples_box_plots_{flank_length}bp.png",
               flank_length=FLANK_LENGTHS),
        # Group1 analysis outputs
        expand(DIVERSITY_DIR / "group1_diversity_metrics_{flank_length}bp.tsv",
               flank_length=FLANK_LENGTHS),
        expand(EVOLUTION_PLOTS_DIR / "frequency_line_plots_{flank_length}bp_present_spectrum.png",
               flank_length=FLANK_LENGTHS),
        expand(EVOLUTION_PLOTS_DIR / "frequency_box_plots_{flank_length}bp_present_spectrum.png",
               flank_length=FLANK_LENGTHS)


rule all_samples_evolution_only:
    """
    Target: All-samples fixation analysis only
    """
    input:
        expand(DIVERSITY_DIR / "all_samples_diversity_metrics_{flank_length}bp.tsv",
               flank_length=FLANK_LENGTHS),
        expand(EVOLUTION_PLOTS_DIR / "all_samples_box_plots_{flank_length}bp.png",
               flank_length=FLANK_LENGTHS)


rule group1_evolution_only:
    """
    Target: Group1 frequency analysis only
    """
    input:
        expand(DIVERSITY_DIR / "group1_diversity_metrics_{flank_length}bp.tsv",
               flank_length=FLANK_LENGTHS),
        expand(EVOLUTION_PLOTS_DIR / "frequency_line_plots_{flank_length}bp_present_spectrum.png",
               flank_length=FLANK_LENGTHS),
        expand(EVOLUTION_PLOTS_DIR / "frequency_box_plots_{flank_length}bp_present_spectrum.png",
               flank_length=FLANK_LENGTHS)


rule consensus_sequences_only:
    """
    Target: Build consensus sequences only (dependency for evolution analysis)
    """
    input:
        expand(
            CONSENSUS_DIR / "{sample}.loci.{side}_flank_{flank_length}bp.consensus.fa",
            sample=ALL_SAMPLES, side=SIDES, flank_length=FLANK_LENGTHS
        )
