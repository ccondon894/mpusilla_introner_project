# ============================================================
# 20_snp_calling.smk - SNP Calling with GATK HaplotypeCaller
# ============================================================
#
# Performs joint SNP calling across all samples using GATK.
# Uses haploid mode (--sample-ploidy 1) for M. pusilla.
#
# Pipeline:
# 1. Index reference genome
# 2. Add read groups to BAMs
# 3. Run HaplotypeCaller per sample (GVCF mode)
# 4. Combine GVCFs
# 5. Joint genotyping
#
# Adapted from:
# - /scratch1/chris/introner-genotyping-final/rules/call_haplotypes_to_CCMP1545.smk
#
# ============================================================

import os
from pathlib import Path

# ============================================================
# CONFIGURATION
# ============================================================

# Output directories
SNP_DIR = RESULTS / "snp_popgen"
VCF_DIR = SNP_DIR / "vcf"
GVCF_DIR = VCF_DIR / "gvcf"
SNP_LOG_DIR = SNP_DIR / "logs"

# Reference genome (using CCMP1545 v3 from config)
SNP_REFERENCE = config["paths"]["references"]["ccmp1545_graph"]

# Source BAM directory (from vg surject alignments)
# These are BAMs aligned to CCMP1545 reference
SOURCE_BAM_DIR = Path("/storage1/chris/introner-genotyping")

# GATK parameters
PLOIDY = config["params"]["gatk"]["ploidy"]

# Samples to call variants on (exclude reference - it's added back as hom-ref in 21_snp_popgen)
SNP_CALL_SAMPLES = [s for s in ALL_SAMPLES if s != REFERENCE]

# ============================================================
# HELPER FUNCTIONS
# ============================================================

def get_source_bam(wildcards):
    """Get source BAM file for a sample"""
    return SOURCE_BAM_DIR / f"{wildcards.sample}.surject_to.CCMP1545.sorted.bam"


# ============================================================
# REFERENCE INDEXING
# ============================================================

rule snp_index_reference:
    """
    Index the reference genome for GATK.

    Creates:
    - FASTA index (.fai)
    - Sequence dictionary (.dict)
    """
    input:
        ref = SNP_REFERENCE
    output:
        fai = VCF_DIR / "reference" / "CCMP1545_v3.fa.fai",
        dict_file = VCF_DIR / "reference" / "CCMP1545_v3.dict",
        ref_copy = VCF_DIR / "reference" / "CCMP1545_v3.fa"
    log:
        SNP_LOG_DIR / "index_reference.log"
    shell:
        """
        mkdir -p {VCF_DIR}/reference
        cp {input.ref} {output.ref_copy}
        samtools faidx {output.ref_copy} 2> {log}
        gatk CreateSequenceDictionary -R {output.ref_copy} -O {output.dict_file} 2>> {log}
        """


# ============================================================
# BAM PREPARATION
# ============================================================

rule add_read_groups:
    """
    Add read group information to BAM files.

    Required by GATK for proper sample identification.
    """
    input:
        bam = get_source_bam
    output:
        bam = VCF_DIR / "bams" / "{sample}.rg.bam"
    params:
        RGID = "{sample}",
        RGLB = "lib1",
        RGPL = "ILLUMINA",
        RGPU = "{sample}.unit1",
        RGSM = "{sample}"
    log:
        SNP_LOG_DIR / "add_read_groups" / "{sample}.log"
    shell:
        """
        mkdir -p {VCF_DIR}/bams
        mkdir -p {SNP_LOG_DIR}/add_read_groups
        gatk AddOrReplaceReadGroups \
            I={input.bam} \
            O={output.bam} \
            RGID={params.RGID} \
            RGLB={params.RGLB} \
            RGPL={params.RGPL} \
            RGPU={params.RGPU} \
            RGSM={params.RGSM} \
            2> {log}
        """


rule index_rg_bam:
    """
    Index BAM files with read groups added.
    """
    input:
        bam = VCF_DIR / "bams" / "{sample}.rg.bam"
    output:
        bai = VCF_DIR / "bams" / "{sample}.rg.bam.bai"
    shell:
        """
        samtools index {input.bam}
        """


# ============================================================
# VARIANT CALLING
# ============================================================

rule haplotype_caller_gvcf:
    """
    Call variants per sample using GATK HaplotypeCaller in GVCF mode.

    Uses haploid mode for M. pusilla.
    GVCF format enables joint genotyping across samples.
    """
    input:
        bam = VCF_DIR / "bams" / "{sample}.rg.bam",
        bai = VCF_DIR / "bams" / "{sample}.rg.bam.bai",
        ref = VCF_DIR / "reference" / "CCMP1545_v3.fa",
        fai = VCF_DIR / "reference" / "CCMP1545_v3.fa.fai",
        dict_file = VCF_DIR / "reference" / "CCMP1545_v3.dict"
    output:
        gvcf = GVCF_DIR / "{sample}.g.vcf.gz",
        gvcf_index = GVCF_DIR / "{sample}.g.vcf.gz.tbi"
    params:
        ploidy = PLOIDY
    threads: 4
    log:
        SNP_LOG_DIR / "haplotype_caller" / "{sample}.log"
    shell:
        """
        mkdir -p {GVCF_DIR}
        mkdir -p {SNP_LOG_DIR}/haplotype_caller
        gatk HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.gvcf} \
            -ERC GVCF \
            --sample-ploidy {params.ploidy} \
            --native-pair-hmm-threads {threads} \
            2> {log}
        """


rule combine_gvcfs:
    """
    Combine individual GVCFs into a single multi-sample GVCF.
    """
    input:
        gvcfs = expand(GVCF_DIR / "{sample}.g.vcf.gz", sample=SNP_CALL_SAMPLES),
        ref = VCF_DIR / "reference" / "CCMP1545_v3.fa"
    output:
        combined = VCF_DIR / "combined.g.vcf.gz"
    params:
        gvcf_args = lambda wildcards, input: " ".join([f"-V {g}" for g in input.gvcfs])
    threads: 4
    log:
        SNP_LOG_DIR / "combine_gvcfs.log"
    shell:
        """
        gatk CombineGVCFs \
            -R {input.ref} \
            {params.gvcf_args} \
            -O {output.combined} \
            2> {log}
        """


rule joint_genotyping:
    """
    Perform joint genotyping on combined GVCF.

    Produces the final joint VCF with all variants.
    """
    input:
        combined = VCF_DIR / "combined.g.vcf.gz",
        ref = VCF_DIR / "reference" / "CCMP1545_v3.fa"
    output:
        vcf = VCF_DIR / "mpusilla.joint.vcf.gz",
        vcf_index = VCF_DIR / "mpusilla.joint.vcf.gz.tbi"
    threads: 8
    log:
        SNP_LOG_DIR / "joint_genotyping.log"
    shell:
        """
        gatk GenotypeGVCFs \
            -R {input.ref} \
            -V {input.combined} \
            -O {output.vcf} \
            2> {log}
        """


rule joint_genotyping_all_sites:
    """
    Perform joint genotyping including non-variant sites.

    Useful for downstream analyses requiring invariant site information.
    """
    input:
        combined = VCF_DIR / "combined.g.vcf.gz",
        ref = VCF_DIR / "reference" / "CCMP1545_v3.fa"
    output:
        vcf = VCF_DIR / "mpusilla.joint.all_sites.vcf.gz",
        vcf_index = VCF_DIR / "mpusilla.joint.all_sites.vcf.gz.tbi"
    threads: 8
    log:
        SNP_LOG_DIR / "joint_genotyping_all_sites.log"
    shell:
        """
        gatk GenotypeGVCFs \
            -R {input.ref} \
            -V {input.combined} \
            -O {output.vcf} \
            --include-non-variant-sites \
            2> {log}
        """


# ============================================================
# TARGET RULES
# ============================================================

rule snp_calling_complete:
    """
    Target: Complete SNP calling pipeline.
    """
    input:
        VCF_DIR / "mpusilla.joint.vcf.gz",
        VCF_DIR / "mpusilla.joint.all_sites.vcf.gz"
