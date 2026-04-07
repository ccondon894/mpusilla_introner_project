# ============================================================
# 11_haplotype_diversity.smk - Haplotype Diversity Ratio Analysis
# ============================================================
#
# Analyzes haplotype diversity ratios (PHDR/AHDR) for introner regions
# using 4-fold degenerate sites as a neutral reference baseline.
#
# Key concepts:
# - PHDR: Present Haplotype Diversity Ratio (π_present / π_4fold)
# - AHDR: Absent Haplotype Diversity Ratio (π_absent / π_4fold)
#
# These ratios compare introner flanking region diversity to
# neutral expectations based on 4-fold degenerate sites,
# controlling for sample set composition.
#
# Steps:
# 1. Parse VCF with 4-fold degenerate sites to create sample-set lookup
# 2. Calculate PHDR/AHDR using sample-set-specific neutral baselines
# 3. Generate visualization plots and summary statistics
#
# Adapted from: /scratch1/chris/introner-genotyping-pipeline/rules/haplotype_diversity_workflow.smk
#
# ============================================================

import os
from pathlib import Path

# ============================================================
# CONFIGURATION
# ============================================================

# Output directories
HAPLOTYPE_DIR = EVOLUTION_DIR / "haplotype_diversity"
HAPLOTYPE_LOOKUP_DIR = HAPLOTYPE_DIR / "lookup_tables"
HAPLOTYPE_RATIOS_DIR = HAPLOTYPE_DIR / "diversity_metrics"
HAPLOTYPE_PLOTS_DIR = HAPLOTYPE_DIR / "plots"
HAPLOTYPE_LOG_DIR = HAPLOTYPE_DIR / "logs"

# Input files
# 4-fold degenerate site VCF (from SNP popgen pipeline or pre-computed)
FOURFOLD_VCF = Path(config["paths"].get("fourfold_vcf",
    "/scratch1/chris/introner-genotyping-pipeline/popgen_snp_data/mpusilla.snps.4d.notMT.group1.allsites.vcf.gz"))

# Introner diversity metrics (from 10_evolution_analysis.smk)
# Uses 200bp flank length by default
INTRONER_METRICS_FILE = DIVERSITY_DIR / "group1_diversity_metrics_200bp.tsv"

# Plot configuration
HAPLOTYPE_FIGSIZE = ",".join(map(str, config["plotting"]["figsize"]["default"]))
HAPLOTYPE_DPI = config["plotting"]["dpi"]


# ============================================================
# RULES
# ============================================================

rule create_4fold_lookup_table:
    """
    Create lookup table for average π values by sample sets from 4-fold VCF.

    Calculates π site-by-site for each unique sample set, then averages
    across all sites. Filters out sites that overlap with introner regions.

    This provides sample-set-specific neutral baselines more accurate than
    pre-computed metrics that assume all samples are present.
    """
    input:
        vcf = FOURFOLD_VCF,
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv"
    output:
        lookup_table = HAPLOTYPE_LOOKUP_DIR / "4fold_avg_pi_present_samples.tsv"
    log:
        HAPLOTYPE_LOG_DIR / "create_4fold_lookup_tables.log"
    shell:
        """
        mkdir -p {HAPLOTYPE_LOOKUP_DIR}
        mkdir -p {HAPLOTYPE_LOG_DIR}

        python {PROJECT_ROOT}/scripts/evolution/create_4fold_vcf_lookup.py \
            --vcf {input.vcf} \
            --genotype_matrix {input.genotype_matrix} \
            --output_dir {HAPLOTYPE_LOOKUP_DIR} \
            2>&1 | tee {log}
        """


rule calculate_haplotype_diversity_ratios:
    """
    Calculate Present and Absent Haplotype Diversity Ratios (PHDR/AHDR).

    For each polymorphic introner:
    - PHDR = π_present_flanks / π_4fold_present_samples
    - AHDR = π_absent_flanks / π_4fold_absent_samples

    Both ratios use sample-set-specific 4D diversity for proper normalization.
    """
    input:
        introner_metrics = INTRONER_METRICS_FILE,
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        fourfold_lookup = HAPLOTYPE_LOOKUP_DIR / "4fold_avg_pi_present_samples.tsv"
    output:
        ratios = HAPLOTYPE_RATIOS_DIR / "haplotype_diversity_ratios.tsv"
    log:
        HAPLOTYPE_LOG_DIR / "calculate_haplotype_diversity_ratios.log"
    shell:
        """
        mkdir -p {HAPLOTYPE_RATIOS_DIR}

        python {PROJECT_ROOT}/scripts/evolution/calculate_haplotype_diversity_ratios.py \
            --introner_metrics {input.introner_metrics} \
            --genotype_matrix {input.genotype_matrix} \
            --fourfold_lookup {input.fourfold_lookup} \
            --output {output.ratios} \
            --min_sites 10 \
            2>&1 | tee {log}
        """


rule plot_haplotype_diversity_ratios:
    """
    Create jitter plots for haplotype diversity deltas with sign tests.

    Generates:
    - Jitter plots of π_introner - π_4d for PHDR/AHDR by frequency
    - Sign test summary (% below neutral baseline at each frequency)
    """
    input:
        diversity_metrics = INTRONER_METRICS_FILE,
        haplotype_ratios = HAPLOTYPE_RATIOS_DIR / "haplotype_diversity_ratios.tsv"
    output:
        plot = HAPLOTYPE_PLOTS_DIR / "haplotype_diversity_ratios_by_frequency.png",
        summary_stats = HAPLOTYPE_RATIOS_DIR / "haplotype_diversity_sign_test_summary.tsv"
    params:
        dpi = HAPLOTYPE_DPI
    log:
        HAPLOTYPE_LOG_DIR / "plot_haplotype_diversity_ratios.log"
    shell:
        """
        mkdir -p {HAPLOTYPE_PLOTS_DIR}

        python {PROJECT_ROOT}/scripts/evolution/plot_haplotype_diversity_ratios.py \
            --diversity_metrics {input.diversity_metrics} \
            --haplotype_ratios {input.haplotype_ratios} \
            --output {output.plot} \
            --summary_stats {output.summary_stats} \
            --dpi {params.dpi} \
            --verbose \
            2>&1 | tee {log}
        """


# ============================================================
# TARGET RULES
# ============================================================

rule all_haplotype_diversity:
    """
    Target: Complete haplotype diversity ratio analysis
    """
    input:
        HAPLOTYPE_PLOTS_DIR / "haplotype_diversity_ratios_by_frequency.png",
        HAPLOTYPE_RATIOS_DIR / "haplotype_diversity_sign_test_summary.tsv"


rule haplotype_lookup_only:
    """
    Target: Generate 4-fold lookup table only
    """
    input:
        HAPLOTYPE_LOOKUP_DIR / "4fold_avg_pi_present_samples.tsv"


rule haplotype_ratios_only:
    """
    Target: Calculate haplotype ratios (without plots)
    """
    input:
        HAPLOTYPE_RATIOS_DIR / "haplotype_diversity_ratios.tsv"
