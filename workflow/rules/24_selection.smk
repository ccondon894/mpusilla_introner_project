# ============================================================
# 24_selection.smk - Selection & Population Structure Analysis
# ============================================================
#
# Performs various population genetics analyses:
# - Tajima's D (neutrality test)
# - Site frequency spectra
# - Principal component analysis (PCA)
# - Linkage disequilibrium (LD)
# - Recombination rate comparisons
#
# Adapted from:
# - /scratch1/chris/introner_vis/basic_popgen/
# - /scratch1/chris/introner_vis/pca/
# - /scratch1/chris/introner_vis/ld/
# - /scratch1/chris/introner_vis/polymorphism_analysis/
# - /scratch1/chris/introner_vis/recombination_analysis/
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
SELECTION_DIR = SNP_DIR / "selection"
PCA_DIR = SELECTION_DIR / "pca"
LD_DIR = SELECTION_DIR / "ld"
SFS_DIR = SELECTION_DIR / "sfs"
RECOMB_DIR = SELECTION_DIR / "recombination"
SELECTION_LOG_DIR = SNP_DIR / "logs" / "selection"

# Pyrho recombination map paths
PYRHO_CCMP1545 = config["paths"]["pyrho"]["ccmp1545"]
PYRHO_RCC1749 = config["paths"]["pyrho"]["rcc1749"]


# ============================================================
# BASIC POPULATION GENETICS STATISTICS
# ============================================================

rule calculate_tajimas_d:
    """
    Calculate Tajima's D statistic for neutrality testing.

    Tajima's D compares two estimates of genetic diversity:
    - Positive D: Balancing selection or population contraction
    - Negative D: Purifying selection or population expansion
    - D ~ 0: Neutral evolution
    """
    input:
        vcf = VCF_DIR / "mpusilla.snps.4d.notMT.group1.vcf.gz"
    output:
        tsv = SELECTION_DIR / "tajimas_d.tsv",
        summary = SELECTION_DIR / "tajimas_d_summary.txt"
    log:
        SELECTION_LOG_DIR / "tajimas_d.log"
    shell:
        """
        mkdir -p {SELECTION_DIR}
        mkdir -p {SELECTION_LOG_DIR}

        python {PROJECT_ROOT}/scripts/popgen/basic_popgen/tajima_d.py \
            --vcf {input.vcf} \
            --output {output.tsv} \
            --summary {output.summary} \
            2> {log}
        """


rule calculate_unfolded_sfs:
    """
    Calculate unfolded site frequency spectrum.

    Uses Group2 (RCC1749, RCC3052) as outgroup to polarize
    ancestral vs derived alleles.
    """
    input:
        vcf = VCF_DIR / "mpusilla.snps.4d.notMT.vcf.gz"
    output:
        sfs = SFS_DIR / "unfolded_1d_sfs.tsv",
        plot = FIGURES_DIR / "unfolded_sfs.pdf"
    params:
        group1_str = ",".join(GROUP1_SAMPLES),
        outgroup_str = ",".join(GROUP2_SAMPLES)
    log:
        SELECTION_LOG_DIR / "unfolded_sfs.log"
    shell:
        """
        mkdir -p {SFS_DIR}

        python {PROJECT_ROOT}/scripts/popgen/basic_popgen/unfolded_1D_afs.py \
            --vcf {input.vcf} \
            --group1 {params.group1_str} \
            --outgroup {params.outgroup_str} \
            --output {output.sfs} \
            --plot {output.plot} \
            2> {log}
        """


rule sfs_density_analysis:
    """
    Perform SFS density analysis comparing introner states.

    Analyzes site frequency spectra in genomic regions with
    vs without introners.
    """
    input:
        vcf = VCF_DIR / "mpusilla.snps.4d.notMT.vcf.gz",
        snpeff_vcf = SNP_DIR / "snpeff" / "mpusilla.snps.snpEff.no_MT.vcf",
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.tsv"
    output:
        tsv = SFS_DIR / "sfs_by_introner_state.tsv",
        plot = FIGURES_DIR / "sfs_density_by_introner.pdf"
    log:
        SELECTION_LOG_DIR / "sfs_density.log"
    shell:
        """
        python {PROJECT_ROOT}/scripts/popgen/polymorphism_analysis/phase2_sfs_analysis.py \
            --vcf {input.vcf} \
            --snpeff {input.snpeff_vcf} \
            --genotype_matrix {input.genotype_matrix} \
            --output {output.tsv} \
            --plot {output.plot} \
            2> {log}
        """


# ============================================================
# PRINCIPAL COMPONENT ANALYSIS
# ============================================================

rule vcf_to_plink:
    """
    Convert VCF to PLINK format for PCA.
    """
    input:
        vcf = VCF_DIR / "mpusilla.snps.4d.notMT.vcf.gz"
    output:
        bed = PCA_DIR / "mpusilla.snps.bed",
        bim = PCA_DIR / "mpusilla.snps.bim",
        fam = PCA_DIR / "mpusilla.snps.fam"
    params:
        prefix = PCA_DIR / "mpusilla.snps"
    log:
        SELECTION_LOG_DIR / "vcf_to_plink.log"
    shell:
        """
        mkdir -p {PCA_DIR}

        plink --vcf {input.vcf} \
            --make-bed \
            --out {params.prefix} \
            --allow-extra-chr \
            --double-id \
            2> {log}
        """


rule ld_pruning:
    """
    LD pruning for PCA to remove correlated SNPs.

    Uses sliding window approach:
    - Window size: 50 SNPs
    - Step size: 5 SNPs
    - r² threshold: 0.2
    """
    input:
        bed = PCA_DIR / "mpusilla.snps.bed",
        bim = PCA_DIR / "mpusilla.snps.bim",
        fam = PCA_DIR / "mpusilla.snps.fam"
    output:
        prune_in = PCA_DIR / "mpusilla.snps.prune.in",
        prune_out = PCA_DIR / "mpusilla.snps.prune.out"
    params:
        prefix = PCA_DIR / "mpusilla.snps",
        window = 50,
        step = 5,
        r2 = 0.2
    log:
        SELECTION_LOG_DIR / "ld_pruning.log"
    shell:
        """
        plink --bfile {params.prefix} \
            --indep-pairwise {params.window} {params.step} {params.r2} \
            --out {params.prefix} \
            --allow-extra-chr \
            2> {log}
        """


rule run_pca:
    """
    Run PCA on LD-pruned SNPs.
    """
    input:
        bed = PCA_DIR / "mpusilla.snps.bed",
        prune_in = PCA_DIR / "mpusilla.snps.prune.in"
    output:
        eigenval = PCA_DIR / "mpusilla.snps.eigenval",
        eigenvec = PCA_DIR / "mpusilla.snps.eigenvec"
    params:
        prefix = PCA_DIR / "mpusilla.snps"
    log:
        SELECTION_LOG_DIR / "pca.log"
    shell:
        """
        plink --bfile {params.prefix} \
            --extract {input.prune_in} \
            --pca 10 \
            --out {params.prefix} \
            --allow-extra-chr \
            2> {log}
        """


rule plot_pca:
    """
    Generate PCA visualization plot.

    Colors samples by group:
    - Group1 (intronerful): Blue
    - Group2 (intronerless): Red
    """
    input:
        eigenval = PCA_DIR / "mpusilla.snps.eigenval",
        eigenvec = PCA_DIR / "mpusilla.snps.eigenvec"
    output:
        pdf = FIGURES_DIR / "pca.pdf",
        png = FIGURES_DIR / "pca.png"
    params:
        group1_str = ",".join(GROUP1_SAMPLES),
        group2_str = ",".join(GROUP2_SAMPLES)
    log:
        SELECTION_LOG_DIR / "plot_pca.log"
    shell:
        """
        python {PROJECT_ROOT}/scripts/popgen/pca/plot_PCA.py \
            --eigenval {input.eigenval} \
            --eigenvec {input.eigenvec} \
            --output {output.pdf} \
            --group1 {params.group1_str} \
            --group2 {params.group2_str} \
            2> {log}

        python {PROJECT_ROOT}/scripts/popgen/pca/plot_PCA.py \
            --eigenval {input.eigenval} \
            --eigenvec {input.eigenvec} \
            --output {output.png} \
            --group1 {params.group1_str} \
            --group2 {params.group2_str} \
            2>> {log}
        """


# ============================================================
# LINKAGE DISEQUILIBRIUM ANALYSIS
# ============================================================

rule compute_ld:
    """
    Compute pairwise linkage disequilibrium (r²) between SNPs.

    Calculates LD within 100kb windows.
    """
    input:
        vcf = VCF_DIR / "mpusilla.snps.4d.notMT.group1.vcf.gz"
    output:
        ld = LD_DIR / "ld_pairwise.tsv.gz"
    params:
        max_distance = 100000  # 100kb
    log:
        SELECTION_LOG_DIR / "compute_ld.log"
    shell:
        """
        mkdir -p {LD_DIR}

        python {PROJECT_ROOT}/scripts/popgen/ld/compute_ld.py \
            --vcf {input.vcf} \
            --output {output.ld} \
            --max_distance {params.max_distance} \
            2> {log}
        """


rule summarize_ld:
    """
    Summarize LD by distance bins.

    Groups SNP pairs into distance bins and calculates
    mean r² per bin for LD decay visualization.
    """
    input:
        ld = LD_DIR / "ld_pairwise.tsv.gz"
    output:
        summary = LD_DIR / "ld_summary_by_distance.tsv"
    params:
        bin_size = 10000  # 10kb bins
    log:
        SELECTION_LOG_DIR / "summarize_ld.log"
    shell:
        """
        python {PROJECT_ROOT}/scripts/popgen/ld/summarize_ld.py \
            {input.ld} \
            {output.summary} \
            --bin_size {params.bin_size} \
            2> {log}
        """


rule plot_ld_decay:
    """
    Plot LD decay with genomic distance.
    """
    input:
        summary = LD_DIR / "ld_summary_by_distance.tsv"
    output:
        pdf = FIGURES_DIR / "ld_decay.pdf",
        png = FIGURES_DIR / "ld_decay.png"
    log:
        SELECTION_LOG_DIR / "plot_ld.log"
    shell:
        """
        python {PROJECT_ROOT}/scripts/popgen/ld/plot_ld.py \
            {input.summary} \
            {output.pdf} \
            2> {log}

        python {PROJECT_ROOT}/scripts/popgen/ld/plot_ld.py \
            {input.summary} \
            {output.png} \
            2>> {log}
        """


# ============================================================
# RECOMBINATION ANALYSIS
# ============================================================

rule compare_recombination_introners:
    """
    Compare recombination rates in introner vs non-introner regions.

    Uses pyrho output and 10kb windows to compare recombination
    landscapes around introner insertion sites.
    """
    input:
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.tsv",
        ccmp1545_gtf = ANNOTATIONS_DIR / "CCMP1545.gtf"
    output:
        tsv = RECOMB_DIR / "recombination_comparison_10kb.tsv",
        summary = RECOMB_DIR / "recombination_comparison_summary.txt"
    params:
        pyrho_dir = PYRHO_CCMP1545,
        window_size = 10000
    log:
        SELECTION_LOG_DIR / "recombination_comparison.log"
    shell:
        """
        mkdir -p {RECOMB_DIR}

        python {PROJECT_ROOT}/scripts/popgen/recombination_analysis/compare_recombination_10kb_windows.py \
            --genotype_matrix {input.genotype_matrix} \
            --pyrho_dir {params.pyrho_dir} \
            --gtf {input.ccmp1545_gtf} \
            --window_size {params.window_size} \
            --output {output.tsv} \
            --summary {output.summary} \
            2> {log}
        """


rule plot_recombination_boxplots:
    """
    Generate boxplots comparing recombination rates by introner category.

    Categories:
    - Non-introner regions
    - All introner regions
    - Polymorphic introners
    - Recent loss introners
    - Recent gain introners
    """
    input:
        tsv = RECOMB_DIR / "recombination_comparison_10kb.tsv"
    output:
        pdf = FIGURES_DIR / "recombination_boxplots.pdf",
        png = FIGURES_DIR / "recombination_boxplots.png"
    log:
        SELECTION_LOG_DIR / "plot_recombination.log"
    shell:
        """
        python {PROJECT_ROOT}/scripts/popgen/recombination_analysis/plot_recombination_boxplots.py \
            --input {input.tsv} \
            --output {output.pdf} \
            2> {log}

        python {PROJECT_ROOT}/scripts/popgen/recombination_analysis/plot_recombination_boxplots.py \
            --input {input.tsv} \
            --output {output.png} \
            2>> {log}
        """


rule recombination_by_frequency:
    """
    Analyze recombination rates by introner frequency category.
    """
    input:
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.tsv"
    output:
        tsv = RECOMB_DIR / "recombination_by_frequency.tsv",
        plot = FIGURES_DIR / "recombination_by_frequency.pdf"
    params:
        pyrho_dir = PYRHO_CCMP1545
    log:
        SELECTION_LOG_DIR / "recombination_by_frequency.log"
    shell:
        """
        python {PROJECT_ROOT}/scripts/popgen/recombination_analysis/compare_recombination_gene_exonic_frequency_based.py \
            --genotype_matrix {input.genotype_matrix} \
            --pyrho_dir {params.pyrho_dir} \
            --output {output.tsv} \
            --plot {output.plot} \
            2> {log}
        """


# ============================================================
# TARGET RULES
# ============================================================

rule selection_complete:
    """
    Target: Complete selection and population structure analysis.
    """
    input:
        # Basic popgen
        SELECTION_DIR / "tajimas_d.tsv",
        SFS_DIR / "unfolded_1d_sfs.tsv",
        # PCA
        PCA_DIR / "mpusilla.snps.eigenvec",
        FIGURES_DIR / "pca.pdf",
        # LD
        LD_DIR / "ld_summary_by_distance.tsv",
        FIGURES_DIR / "ld_decay.pdf",
        # Recombination
        RECOMB_DIR / "recombination_comparison_10kb.tsv",
        FIGURES_DIR / "recombination_boxplots.pdf"


rule basic_popgen_only:
    """
    Target: Basic population genetics statistics only.
    """
    input:
        SELECTION_DIR / "tajimas_d.tsv",
        SFS_DIR / "unfolded_1d_sfs.tsv",
        FIGURES_DIR / "unfolded_sfs.pdf"


rule pca_only:
    """
    Target: PCA analysis only.
    """
    input:
        PCA_DIR / "mpusilla.snps.eigenvec",
        FIGURES_DIR / "pca.pdf"


rule ld_only:
    """
    Target: LD analysis only.
    """
    input:
        LD_DIR / "ld_summary_by_distance.tsv",
        FIGURES_DIR / "ld_decay.pdf"


rule recombination_only:
    """
    Target: Recombination analysis only.
    """
    input:
        RECOMB_DIR / "recombination_comparison_10kb.tsv",
        FIGURES_DIR / "recombination_boxplots.pdf"
