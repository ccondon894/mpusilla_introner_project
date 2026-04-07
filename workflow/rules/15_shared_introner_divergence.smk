# ============================================================================
# 15_shared_introner_divergence.smk
# Test whether shared introners (present in both Group 1 and Group 2) are
# ancestral or convergent insertions by comparing introner body Dxy to
# flanking region Dxy.
# ============================================================================

import os
from pathlib import Path

# Output directories
SHARED_DIV_DIR = EVOLUTION_DIR / "shared_introner_divergence"
SHARED_DIV_CONSENSUS = SHARED_DIV_DIR / "consensus_sequences"
SHARED_DIV_ALIGNMENTS = SHARED_DIV_DIR / "alignments"
SHARED_DIV_METRICS = SHARED_DIV_DIR / "metrics"
SHARED_DIV_FIGURES = FIGURES_DIR / "shared_introner_divergence"

# Script directory
SHARED_DIV_SCRIPTS = str(Path(workflow.basedir).parent / "scripts" / "shared_introner_divergence")

# Input directories (reuse from rule 10)
COVERAGE_BAM_DIR = GENOTYPING_DIR / "coverage" / "bams"

# Non-reference samples (need mpileup consensus pipeline)
NON_REF_SAMPLES = [s for s in ALL_SAMPLES if s != REFERENCE]

# Region types for shared introner analysis
SHARED_REGIONS = ["introner_body", "left_flank", "right_flank"]


# ============================================================================
# TARGET RULE
# ============================================================================

rule shared_introner_divergence_all:
    """
    Target: Compare introner body divergence to flanking divergence
    for shared loci to distinguish ancestral from convergent insertions.
    """
    input:
        SHARED_DIV_METRICS / "introner_vs_flank_dxy_all_shared.tsv",
        SHARED_DIV_METRICS / "introner_vs_flank_dxy_fixed_shared.tsv",
        SHARED_DIV_METRICS / "introner_vs_flank_dxy_polymorphic_shared.tsv",
        SHARED_DIV_FIGURES / "dxy_scatter_all_shared.png"


# ============================================================================
# STEP 1: IDENTIFY SHARED LOCI AND GENERATE BED FILES
# ============================================================================

rule identify_shared_loci:
    """Identify shared introner loci and generate BED files for body + flanks."""
    input:
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv"
    output:
        classification = SHARED_DIV_DIR / "shared_loci_classification.json",
        beds = expand(
            SHARED_DIV_CONSENSUS / "{sample}.shared.{region}.bed",
            sample=ALL_SAMPLES,
            region=SHARED_REGIONS
        )
    params:
        group1 = ' '.join(GROUP1_SAMPLES),
        group2 = ' '.join(GROUP2_SAMPLES),
        output_dir = str(SHARED_DIV_CONSENSUS)
    log:
        SHARED_DIV_DIR / "logs" / "identify_shared_loci.log"
    shell:
        """
        mkdir -p $(dirname {log}) {params.output_dir} && \
        python3 {SHARED_DIV_SCRIPTS}/identify_shared_loci.py \
            --genotype-matrix {input.genotype_matrix} \
            --group1-samples {params.group1} \
            --group2-samples {params.group2} \
            --flank-length 100 \
            --output-dir {params.output_dir} \
            --classification-json {output.classification} \
            2>&1 | tee {log}
        """


# ============================================================================
# STEP 2: BUILD CONSENSUS SEQUENCES (mpileup pipeline)
# ============================================================================

rule shared_introner_mpileup:
    """Generate mpileup for shared introner regions (non-reference samples)."""
    input:
        bam = COVERAGE_BAM_DIR / "{sample}.sorted.bam",
        bed = SHARED_DIV_CONSENSUS / "{sample}.shared.{region}.bed",
        fa = ASSEMBLIES_DIR / "{sample}.vg_paths.fa",
    output:
        mpileup = temp(SHARED_DIV_CONSENSUS / "{sample}.shared.{region}.mpileup.bcf")
    params:
        temp_dir = lambda wc: str(SHARED_DIV_CONSENSUS / f"tmp_{wc.sample}_{wc.region}_sort")
    wildcard_constraints:
        region = "introner_body|left_flank|right_flank",
        sample = "|".join(NON_REF_SAMPLES)
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


rule shared_introner_call_variants:
    """Call variants from mpileup for shared introner consensus building."""
    input:
        mpileup = SHARED_DIV_CONSENSUS / "{sample}.shared.{region}.mpileup.bcf",
        fa = ASSEMBLIES_DIR / "{sample}.vg_paths.fa"
    output:
        vcf = temp(SHARED_DIV_CONSENSUS / "{sample}.shared.{region}.vcf.gz"),
        filt_vcf = temp(SHARED_DIV_CONSENSUS / "{sample}.shared.{region}.filtered.vcf.gz")
    wildcard_constraints:
        region = "introner_body|left_flank|right_flank",
        sample = "|".join(NON_REF_SAMPLES)
    shell:
        """
        bcftools call -m --ploidy 1 -Ou {input.mpileup} | \
        bcftools norm -f {input.fa} -m +both -Oz -o {output.vcf}

        bcftools filter -i 'DP>=10' {output.vcf} -o {output.filt_vcf}
        tabix -p vcf {output.filt_vcf}
        """


rule shared_introner_build_consensus:
    """Build consensus sequences for shared introner regions (non-reference samples)."""
    input:
        filt_vcf = SHARED_DIV_CONSENSUS / "{sample}.shared.{region}.filtered.vcf.gz",
        bed = SHARED_DIV_CONSENSUS / "{sample}.shared.{region}.bed",
        fa = ASSEMBLIES_DIR / "{sample}.vg_paths.fa"
    output:
        cons = SHARED_DIV_CONSENSUS / "{sample}.shared.{region}.consensus.fa"
    params:
        cons_genome = lambda wc: str(SHARED_DIV_CONSENSUS / f"{wc.sample}.shared.{wc.region}.tmp_consensus.fa"),
        loci = lambda wc: str(SHARED_DIV_CONSENSUS / f"{wc.sample}.shared.{wc.region}.tmp_loci.txt")
    wildcard_constraints:
        region = "introner_body|left_flank|right_flank",
        sample = "|".join(NON_REF_SAMPLES)
    shell:
        """
        bcftools consensus -a N -f {input.fa} -o {params.cons_genome} {input.filt_vcf}
        samtools faidx {params.cons_genome}
        awk '{{print $1":"$2"-"$3}}' {input.bed} > {params.loci}
        samtools faidx {params.cons_genome} $(cat {params.loci}) > {output.cons}
        rm -f {params.cons_genome} {params.cons_genome}.fai {params.loci}
        """


rule shared_introner_ref_consensus:
    """Extract reference consensus sequences for shared introner regions (no variant calling)."""
    input:
        bed = SHARED_DIV_CONSENSUS / f"{REFERENCE}.shared.{{region}}.bed",
        fa = ASSEMBLIES_DIR / f"{REFERENCE}.vg_paths.fa",
    output:
        fa = SHARED_DIV_CONSENSUS / f"{REFERENCE}.shared.{{region}}.consensus.fa"
    wildcard_constraints:
        region = "introner_body|left_flank|right_flank"
    shell:
        """
        bedtools getfasta -fi {input.fa} -bed {input.bed} -fo {output.fa}
        """


# ============================================================================
# STEP 3: BUILD MAFFT ALIGNMENTS
# ============================================================================

rule align_shared_introners:
    """Build MAFFT alignments for all shared introner loci."""
    input:
        classification = SHARED_DIV_DIR / "shared_loci_classification.json",
        # Require all consensus files
        ref_consensus = expand(
            SHARED_DIV_CONSENSUS / f"{REFERENCE}.shared.{{region}}.consensus.fa",
            region=SHARED_REGIONS
        ),
        sample_consensus = expand(
            SHARED_DIV_CONSENSUS / "{sample}.shared.{region}.consensus.fa",
            sample=NON_REF_SAMPLES,
            region=SHARED_REGIONS
        )
    output:
        sentinel = touch(SHARED_DIV_ALIGNMENTS / ".alignments_complete")
    params:
        consensus_dir = str(SHARED_DIV_CONSENSUS),
        output_dir = str(SHARED_DIV_ALIGNMENTS)
    log:
        SHARED_DIV_DIR / "logs" / "align_shared_introners.log"
    shell:
        """
        mkdir -p {params.output_dir} $(dirname {log}) && \
        python3 {SHARED_DIV_SCRIPTS}/build_shared_introner_alignments.py \
            --classification-json {input.classification} \
            --consensus-dir {params.consensus_dir} \
            --output-dir {params.output_dir} \
            2>&1 | tee {log}
        """


# ============================================================================
# STEP 4: CALCULATE DXY AND COMPARE
# ============================================================================

rule calculate_introner_flank_dxy:
    """Calculate and compare introner body Dxy vs flanking Dxy."""
    input:
        classification = SHARED_DIV_DIR / "shared_loci_classification.json",
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        alignment_sentinel = SHARED_DIV_ALIGNMENTS / ".alignments_complete"
    output:
        all_tsv = SHARED_DIV_METRICS / "introner_vs_flank_dxy_all_shared.tsv",
        fixed_tsv = SHARED_DIV_METRICS / "introner_vs_flank_dxy_fixed_shared.tsv",
        poly_tsv = SHARED_DIV_METRICS / "introner_vs_flank_dxy_polymorphic_shared.tsv",
        summary = SHARED_DIV_METRICS / "shared_introner_divergence_summary.txt",
        scatter = SHARED_DIV_FIGURES / "dxy_scatter_all_shared.png"
    params:
        alignment_dir = str(SHARED_DIV_ALIGNMENTS),
        output_dir = str(SHARED_DIV_METRICS),
        figures_dir = str(SHARED_DIV_FIGURES),
        group1 = ' '.join(GROUP1_SAMPLES),
        group2 = ' '.join(GROUP2_SAMPLES)
    log:
        SHARED_DIV_DIR / "logs" / "calculate_introner_flank_dxy.log"
    shell:
        """
        mkdir -p {params.output_dir} {params.figures_dir} $(dirname {log}) && \
        python3 {SHARED_DIV_SCRIPTS}/calculate_introner_vs_flank_dxy.py \
            --alignment-dir {params.alignment_dir} \
            --classification-json {input.classification} \
            --genotype-matrix {input.genotype_matrix} \
            --group1-samples {params.group1} \
            --group2-samples {params.group2} \
            --output-dir {params.output_dir} \
            --figures-dir {params.figures_dir} \
            2>&1 | tee {log}
        """
