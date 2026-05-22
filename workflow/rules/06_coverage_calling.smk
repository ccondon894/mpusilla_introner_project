# ============================================================
# 06_coverage_calling.smk - Coverage-Based Introner Calling
# ============================================================
#
# Validates and refines introner genotype calls using read coverage.
# This step improves confidence in presence/absence calls by
# examining read depth across introner loci.
#
# Steps:
# 1. Align reads to sample assemblies (if not already done)
# 2. Create BED files of introner loci from genotype matrix
# 3. Calculate read depth at each locus (samtools depth)
# 4. Call introner presence/absence based on coverage patterns
# 5. Update genotype matrix with coverage-validated calls
#
# Adapted from: /scratch1/chris/introner-genotyping-pipeline/rules/call_introner_presence.smk
#
# ============================================================

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
    """
    Align sample reads to its own graph assembly.

    Uses BWA-MEM for alignment, producing a coordinate-sorted BAM file.
    This is used to calculate read depth at introner loci.
    """
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
    """
    Create BED file of introner loci for a specific sample.

    Extracts loci coordinates from the genotype matrix for use
    with samtools depth calculation.
    """
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
    """
    Calculate read depth at each introner locus using samtools depth.

    Options:
    -Q: Minimum mapping quality (default 30)
    -a: Output all positions (including zero depth)
    -J: Include reads with deletions in depth calculation
    -b: BED file specifying regions
    """
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
    """
    Call introner presence/absence based on coverage patterns.

    Examines read depth across each locus to determine:
    - Present (1): Consistent coverage across entire locus
    - Absent (2): Coverage gap where introner would be
    - Missing data (3): Insufficient coverage for confident call

    The script uses coverage thresholds to make these determinations.
    """
    input:
        bed = COVERAGE_DIR / "{sample}.loci.bed",
        depth = COVERAGE_DIR / "{sample}.locus_depth.tsv",
        gt_matrix = GENOTYPING_DIR / "genotype_matrix.tsv"
    output:
        calls = COVERAGE_DIR / "{sample}.loci.coverage_calls.bed"
    # log:
    #     COVERAGE_LOG_DIR / "{sample}.coverage_calls.log"
    shell:
        # """
        # python {PROJECT_ROOT}/scripts/genotyping/call_locus_from_bam_coverage_v2.py \
        #     {input.bed} {input.depth} {input.gt_matrix} {output.calls} > {log} 2>&1
        # """
        """
        introner-caller {input.bed} {input.depth} {output.calls}
        """


rule update_genotype_matrix_with_coverage:
    """
    Update genotype matrix with coverage-validated calls.

    Incorporates coverage-based calls into the genotype matrix,
    potentially upgrading missing data (3) to presence (1) or
    absence (2) calls based on read evidence.

    Also appends derived, paper-facing classification columns that separate:
      - per-group counts and callability
      - per-group presence pattern
      - within-group orthology confidence
      - cross-group origin and confidence

    The legacy within_group_status / cross_group_status columns are retained
    for downstream compatibility.
    """
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
    """
    Generate summary statistics for the final genotype matrix.
    """
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
    """
    Plot introner family distributions across all samples.
    """
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
    """
    Plot the introner-present 2D AFS from the final genotype matrix.

    Excludes independent insertion classes from the cross-group comparison.
    """
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
    """
    Target: Complete coverage-based validation
    """
    input:
        GENOTYPING_DIR / "genotype_matrix.final.tsv",
        GENOTYPING_DIR / "genotype_matrix_summary.txt",
        FIGURES_DIR / "genotyping" / "introner_family_distributions.png",
        FIGURES_DIR / "genotyping" / "introner_2d_afs_present_count.png"


rule coverage_calls_only:
    """
    Target: Generate coverage calls without updating matrix
    """
    input:
        expand(COVERAGE_DIR / "{sample}.loci.coverage_calls.bed",
               sample=NON_REFERENCE_SAMPLES)


rule bams_only:
    """
    Target: Generate all BAM alignments only
    """
    input:
        expand(BAM_DIR / "{sample}.sorted.bam",
               sample=NON_REFERENCE_SAMPLES)
