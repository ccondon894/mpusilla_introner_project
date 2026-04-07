# ============================================================
# 21_snp_popgen.smk - SNP Filtering and 4-Fold Degenerate Sites
# ============================================================
#
# Filters the joint VCF and extracts 4-fold degenerate sites
# for population genetics analysis.
#
# Pipeline:
# 1. Run degenotate to identify 4-fold degenerate sites
# 2. Extract SNPs from joint VCF
# 3. Add reference sample (CCMP1545)
# 4. Filter for 4-fold degenerate sites
# 5. Create sample subsets (intronerful, no missing data)
# 6. Run SnpEff annotation
#
# Adapted from:
# - /scratch1/chris/introner_vis/data_prep/commands.sh
# - /scratch1/chris/introner_vis/degenotate/degenotate.sh
#
# ============================================================

import os
from pathlib import Path

# ============================================================
# CONFIGURATION
# ============================================================

# Input/output directories
SNP_DIR = RESULTS / "snp_popgen"
VCF_DIR = SNP_DIR / "vcf"
DEGENOTATE_DIR = SNP_DIR / "degenotate"
SNPEFF_DIR = SNP_DIR / "snpeff"
SNP_LOG_DIR = SNP_DIR / "logs"

# Reference files
SNP_REFERENCE = config["paths"]["references"]["ccmp1545_graph"]
REFERENCE_GTF = config["paths"]["references"]["ccmp1545_gtf"]

# Excluded scaffold (mitochondrial/problematic)
EXCLUDED_SCAFFOLD = config["mating_type_region"]["scaffold"]  # scaffold_2

# Sample lists for filtering
INTRONERFUL_SAMPLES = GROUP1_SAMPLES  # Group1 samples have introners
GROUP2_SAMPLES_LIST = GROUP2_SAMPLES  # RCC1749, RCC3052


# ============================================================
# DEGENOTATE - 4-FOLD DEGENERATE SITE IDENTIFICATION
# ============================================================

rule run_degenotate:
    """
    Run degenotate to identify 4-fold degenerate sites.

    Uses the reference genome and GTF annotation to find
    synonymous codon positions.
    """
    input:
        genome = SNP_REFERENCE,
        gtf = ANNOTATIONS_DIR / "CCMP1545.gtf"
    output:
        bed = DEGENOTATE_DIR / "degeneracy-all-sites.bed",
        bed_4d = DEGENOTATE_DIR / "degeneracy-all-sites.4d.bed.gz"
    params:
        outdir = DEGENOTATE_DIR,
        degeneracy_level = 4
    log:
        SNP_LOG_DIR / "degenotate.log"
    shell:
        """
        mkdir -p {params.outdir}

        # Run degenotate
        degenotate.py \
            -g {input.genome} \
            -x {params.degeneracy_level} \
            -a {input.gtf} \
            -o {params.outdir} \
            --overwrite \
            2> {log}

        # Extract 4-fold degenerate sites and compress
        grep -P "\\t4\\t" {output.bed} | gzip > {output.bed_4d}
        """


# ============================================================
# VCF FILTERING AND PROCESSING
# ============================================================

rule extract_snps:
    """
    Extract SNPs from joint VCF (exclude indels and multi-allelic sites).

    Filters for:
    - Single nucleotide changes only (REF and ALT both length 1)
    """
    input:
        vcf = VCF_DIR / "mpusilla.joint.vcf.gz"
    output:
        vcf = VCF_DIR / "mpusilla.snps.vcf.gz"
    shell:
        """
        (
            gzip -dc {input.vcf} | grep -e '^#'
            gzip -dc {input.vcf} | grep -v '^#' | awk 'length($4) == 1 && length($5) == 1'
        ) | gzip > {output.vcf}
        """


rule add_reference_sample:
    """
    Add CCMP1545 as a sample column with homozygous reference genotypes.

    The reference genome (CCMP1545) isn't in the VCF as a sample since
    reads were aligned to it. This adds it for complete analysis.
    """
    input:
        vcf = VCF_DIR / "mpusilla.snps.vcf.gz",
        awk_script = PROJECT_ROOT / "scripts" / "popgen" / "insert_reference.awk"
    output:
        vcf = VCF_DIR / "mpusilla.snps.with_ref.vcf.gz"
    shell:
        """
        gzip -dc {input.vcf} | {input.awk_script} | gzip > {output.vcf}
        """


rule extract_4fold_sites:
    """
    Filter VCF for 4-fold degenerate sites only.

    Also excludes scaffold_2 (mating type region / MT).
    """
    input:
        vcf = VCF_DIR / "mpusilla.snps.with_ref.vcf.gz",
        bed_4d = DEGENOTATE_DIR / "degeneracy-all-sites.4d.bed.gz",
        script = PROJECT_ROOT / "scripts" / "popgen" / "extract_4d.py"
    output:
        vcf = VCF_DIR / "mpusilla.snps.4d.notMT.vcf.gz"
    params:
        exclude_scaffold = EXCLUDED_SCAFFOLD
    shell:
        """
        python {input.script} \
            {input.bed_4d} \
            {input.vcf} \
            | grep -v '{params.exclude_scaffold}' \
            | bgzip > {output.vcf}
        """


rule create_intronerful_subset:
    """
    Create subset VCF with only intronerful (Group1) samples.

    Extracts columns for Group1 samples plus CCMP1545 reference.
    """
    input:
        vcf = VCF_DIR / "mpusilla.snps.4d.notMT.vcf.gz"
    output:
        vcf = VCF_DIR / "mpusilla.snps.4d.notMT.intronerful.vcf.gz"
    params:
        # Build column selection string for Group1 samples
        # First 9 columns are VCF format, then samples
        samples = ",".join(GROUP1_SAMPLES)
    shell:
        """
        # Get header
        gzip -dc {input.vcf} | grep '^##' > {output.vcf}.tmp

        # Get column header and filter for samples
        gzip -dc {input.vcf} | grep '^#CHROM' | head -1 > {output.vcf}.header

        # Use bcftools to subset samples
        bcftools view -s {params.samples} {input.vcf} | bgzip > {output.vcf}
        rm -f {output.vcf}.tmp {output.vcf}.header
        """


rule filter_no_missing_data:
    """
    Create subset VCF with no missing genotypes.

    Filters for sites where all samples have calls (AN = 2*n_samples for diploid,
    or AN = n_samples for haploid).
    """
    input:
        vcf = VCF_DIR / "mpusilla.snps.4d.notMT.vcf.gz"
    output:
        vcf = VCF_DIR / "mpusilla.snps.4d.notMT.nomissing.vcf.gz"
    params:
        # For haploid, AN should equal number of samples
        expected_an = len(ALL_SAMPLES)
    shell:
        """
        (
            gzip -dc {input.vcf} | grep -e '^#'
            gzip -dc {input.vcf} | grep -v '^#' | grep 'AN={params.expected_an}'
        ) | gzip > {output.vcf}
        """


rule create_group1_subset:
    """
    Create subset VCF with only Group1 samples for within-group analysis.
    """
    input:
        vcf = VCF_DIR / "mpusilla.snps.4d.notMT.vcf.gz"
    output:
        vcf = VCF_DIR / "mpusilla.snps.4d.notMT.group1.vcf.gz"
    params:
        samples = ",".join(GROUP1_SAMPLES)
    shell:
        """
        bcftools view -s {params.samples} {input.vcf} | bgzip > {output.vcf}
        """


# ============================================================
# SNPEFF ANNOTATION
# ============================================================

rule build_snpeff_database:
    """
    Build SnpEff database for CCMP1545.

    Creates a custom database from the reference genome and GTF.
    """
    input:
        genome = SNP_REFERENCE,
        gtf = ANNOTATIONS_DIR / "CCMP1545.gtf"
    output:
        done = SNPEFF_DIR / "snpeff_db.done"
    params:
        db_name = "CCMP1545_v3",
        data_dir = SNPEFF_DIR / "data"
    log:
        SNP_LOG_DIR / "snpeff_build.log"
    shell:
        """
        mkdir -p {params.data_dir}/{params.db_name}

        # Copy files with required names
        cp {input.genome} {params.data_dir}/{params.db_name}/sequences.fa
        cp {input.gtf} {params.data_dir}/{params.db_name}/genes.gtf

        # Build database
        snpEff build -gtf22 -v {params.db_name} -dataDir {params.data_dir} 2> {log}

        touch {output.done}
        """


rule run_snpeff:
    """
    Annotate SNPs with functional effects using SnpEff.
    """
    input:
        vcf = VCF_DIR / "mpusilla.snps.vcf.gz",
        db_done = SNPEFF_DIR / "snpeff_db.done"
    output:
        vcf = SNPEFF_DIR / "mpusilla.snps.snpEff.vcf",
        html = SNPEFF_DIR / "snpEff_summary.html",
        genes = SNPEFF_DIR / "snpEff_genes.txt"
    params:
        db_name = "CCMP1545_v3",
        data_dir = SNPEFF_DIR / "data"
    log:
        SNP_LOG_DIR / "snpeff_annotate.log"
    shell:
        """
        gzip -dc {input.vcf} | \
        snpEff ann \
            -v {params.db_name} \
            -dataDir {params.data_dir} \
            -stats {output.html} \
            > {output.vcf} \
            2> {log}

        mv snpEff_genes.txt {output.genes} 2>/dev/null || true
        """


rule summarize_site_classes:
    """
    Summarize SNP counts by site class (non-synonymous, synonymous, 4-fold, non-CDS).

    Excludes scaffold_2 (mating type region).
    """
    input:
        bed = DEGENOTATE_DIR / "degeneracy-all-sites.bed",
        vcf = VCF_DIR / "mpusilla.snps.with_ref.vcf.gz",
        script = PROJECT_ROOT / "scripts" / "popgen" / "basic_popgen" / "site_class_summary.py"
    output:
        tsv = SNP_DIR / "site_class_summary.tsv"
    params:
        exclude_scaffold = EXCLUDED_SCAFFOLD
    log:
        SNP_LOG_DIR / "site_class_summary.log"
    shell:
        """
        python {input.script} \
            --bed {input.bed} \
            --vcf {input.vcf} \
            --output {output.tsv} \
            --exclude-scaffold {params.exclude_scaffold} \
            2> {log}
        """


rule filter_snpeff_no_mt:
    """
    Remove mating type region from SnpEff annotated VCF.
    """
    input:
        vcf = SNPEFF_DIR / "mpusilla.snps.snpEff.vcf"
    output:
        vcf = SNPEFF_DIR / "mpusilla.snps.snpEff.no_MT.vcf"
    params:
        exclude_scaffold = EXCLUDED_SCAFFOLD
    shell:
        """
        grep -v '{params.exclude_scaffold}' {input.vcf} > {output.vcf}
        """


# ============================================================
# TARGET RULES
# ============================================================

rule snp_popgen_complete:
    """
    Target: Complete SNP filtering and annotation.
    """
    input:
        VCF_DIR / "mpusilla.snps.4d.notMT.vcf.gz",
        VCF_DIR / "mpusilla.snps.4d.notMT.intronerful.vcf.gz",
        VCF_DIR / "mpusilla.snps.4d.notMT.nomissing.vcf.gz",
        VCF_DIR / "mpusilla.snps.4d.notMT.group1.vcf.gz",
        SNPEFF_DIR / "mpusilla.snps.snpEff.vcf",
        SNPEFF_DIR / "mpusilla.snps.snpEff.no_MT.vcf"


rule degenotate_only:
    """
    Target: Run degenotate only.
    """
    input:
        DEGENOTATE_DIR / "degeneracy-all-sites.4d.bed.gz"


rule fourD_vcf_only:
    """
    Target: Generate 4-fold degenerate VCF only.
    """
    input:
        VCF_DIR / "mpusilla.snps.4d.notMT.vcf.gz"
