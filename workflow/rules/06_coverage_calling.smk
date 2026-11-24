# 06_coverage_calling.smk - Coverage-Based Introner Calling

import os
from pathlib import Path

# Config

# Output directories
COVERAGE_DIR = GENOTYPING_DIR / "coverage"
COVERAGE_LOG_DIR = COVERAGE_DIR / "logs"
BAM_DIR = COVERAGE_DIR / "bams"

# Input directories
READS_DIR = Path(config["paths"]["original"]["reads"])

# Parameters from config
MAPQ_THRESHOLD = config["params"]["coverage"]["mapq_threshold"]
MIN_DEPTH = config["params"]["coverage"]["min_depth"]
BWA_THREADS = config["params"]["bwa"]["threads"]

# Non-reference samples (we align reads for these)
# CCMP1545 is the reference so we skip it for coverage analysis
NON_REFERENCE_SAMPLES = [s for s in ALL_SAMPLES if s != REFERENCE]

# Rules
rule align_reads_to_assembly:
    """Align sample reads to its own graph assembly."""
    input:
        fa = ASSEMBLIES_DIR / "{sample}.vg_paths.fa",
        fq1 = lambda wildcards: READS_DIR / wildcards.sample / f"{wildcards.sample}.R1.sorted.fastq",
        fq2 = lambda wildcards: READS_DIR / wildcards.sample / f"{wildcards.sample}.R2.sorted.fastq"
    output:
        bam = BAM_DIR / "{sample}.sorted.bam",
        bai = BAM_DIR / "{sample}.sorted.bam.bai"
    log:
        COVERAGE_LOG_DIR / "{sample}.align_reads.log"
    threads: BWA_THREADS
    shell:
        """
        mkdir -p {BAM_DIR}
        mkdir -p {COVERAGE_LOG_DIR}

        bwa mem -t {threads} {input.fa} {input.fq1} {input.fq2} 2> {log} | \
        samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """


rule make_loci_beds:
    """Create BED file of introner loci for a specific sample."""
    input:
        gt_matrix = GENOTYPING_DIR / "genotype_matrix.tsv"
    output:
        bed = COVERAGE_DIR / "{sample}.loci.bed"
    shell:
        """
        mkdir -p {COVERAGE_DIR}
        python {PROJECT_ROOT}/scripts/genotyping/make_sample_bed_from_gt_matrix.py \
            {input.gt_matrix} {wildcards.sample} {output.bed}
        """


rule calculate_locus_depth:
    """Calculate read depth at each introner locus using samtools depth."""
    input:
        bed = COVERAGE_DIR / "{sample}.loci.bed",
        bam = BAM_DIR / "{sample}.sorted.bam"
    output:
        depth = COVERAGE_DIR / "{sample}.locus_depth.tsv"
    params:
        mapq = MAPQ_THRESHOLD
    shell:
        """
        # Skip header/track lines in BED file
        awk '!/^#/ && !/^track/ && !/^browser/ && !/^contig/' {input.bed} | \
        samtools depth -Q {params.mapq} -a -J -b - {input.bam} > {output.depth}
        """


rule call_introner_presence:
    """Call introner presence/absence based on coverage patterns."""
    input:
        bed = COVERAGE_DIR / "{sample}.loci.bed",
        depth = COVERAGE_DIR / "{sample}.locus_depth.tsv",
        gt_matrix = GENOTYPING_DIR / "genotype_matrix.tsv"
    output:
        calls = COVERAGE_DIR / "{sample}.loci.coverage_calls.bed"
    shell:
        """
        introner-caller {input.bed} {input.depth} {output.calls}
        """


rule update_genotype_matrix_with_coverage:
    """Update genotype matrix with coverage-validated calls."""
    input:
        gt_matrix = GENOTYPING_DIR / "genotype_matrix.tsv",
        coverage_calls = expand(COVERAGE_DIR / "{sample}.loci.coverage_calls.bed",
                                sample=NON_REFERENCE_SAMPLES)
    output:
        updated_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv"
    params:
        coverage_dir = COVERAGE_DIR
    shell:
        """
        python {PROJECT_ROOT}/scripts/genotyping/update_gt_matrixes_bam_cov.py \
            {input.gt_matrix} {params.coverage_dir} {output.updated_matrix}
        """


rule genotype_matrix_summary:
    """Generate summary statistics for the final genotype matrix."""
    input:
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv"
    output:
        summary = GENOTYPING_DIR / "genotype_matrix_summary.txt"
    params:
        group1 = ' '.join(GROUP1_SAMPLES),
        group2 = ' '.join(GROUP2_SAMPLES)
    shell:
        """
        python {PROJECT_ROOT}/scripts/genotyping/genotype_matrix_summary.py \
            --genotype-matrix {input.matrix} \
            --group1-samples {params.group1} \
            --group2-samples {params.group2} \
            --output {output.summary}
        """


rule introner_family_distribution_plot:
    """Plot introner family distributions across all samples."""
    input:
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv"
    output:
        plot = FIGURES_DIR / "genotyping" / "introner_family_distributions.png"
    params:
        group1 = ' '.join(GROUP1_SAMPLES),
        group2 = ' '.join(GROUP2_SAMPLES)
    shell:
        """
        mkdir -p $(dirname {output.plot})
        python {PROJECT_ROOT}/scripts/genotyping/introner_family_distributions_panel.py \
            --genotype-matrix {input.matrix} \
            --group1-samples {params.group1} \
            --group2-samples {params.group2} \
            --output {output.plot} \
            --proportions
        """


rule introner_2d_afs_plot:
    """Plot the introner-present 2D AFS from the final genotype matrix."""
    input:
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv"
    output:
        pdf = FIGURES_DIR / "genotyping" / "introner_2d_afs_present_count.pdf",
        png = FIGURES_DIR / "genotyping" / "introner_2d_afs_present_count.png",
        tsv = GENOTYPING_DIR / "introner_2d_afs_present_count.tsv",
        summary = GENOTYPING_DIR / "introner_2d_afs_present_count_summary.txt"
    params:
        group1 = ' '.join(GROUP1_SAMPLES),
        group2 = ' '.join(GROUP2_SAMPLES)
    shell:
        """
        mkdir -p $(dirname {output.png}) $(dirname {output.tsv})
        python {PROJECT_ROOT}/scripts/genotyping/plot_introner_2d_afs.py \
            --genotype-matrix {input.matrix} \
            --group1-samples {params.group1} \
            --group2-samples {params.group2} \
            --output-pdf {output.pdf} \
            --output-png {output.png} \
            --output-tsv {output.tsv} \
            --summary {output.summary}
        """


# Target Rules

rule all_coverage_calling:
    """Target: Complete coverage-based validation"""
    input:
        GENOTYPING_DIR / "genotype_matrix.final.tsv",
        GENOTYPING_DIR / "genotype_matrix_summary.txt",
        FIGURES_DIR / "genotyping" / "introner_family_distributions.png",
        FIGURES_DIR / "genotyping" / "introner_2d_afs_present_count.png"


rule coverage_calls_only:
    """Target: Generate coverage calls without updating matrix"""
    input:
        expand(COVERAGE_DIR / "{sample}.loci.coverage_calls.bed",
               sample=NON_REFERENCE_SAMPLES)


rule bams_only:
    """Target: Generate all BAM alignments only"""
    input:
        expand(BAM_DIR / "{sample}.sorted.bam",
               sample=NON_REFERENCE_SAMPLES)
