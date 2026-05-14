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
SNPEFF_DIR = SNP_DIR / "snpeff"

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
    conda: "../envs/popgen.yaml"
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


rule plot_sfs_by_class:
    """
    Plots the SFS frequencies for synonymous/nonsynonymous and introner allele frequencies
    """
    input:
        snpeff_vcf = SNPEFF_DIR / "mpusilla.snps.snpEff.no_MT.vcf",
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        script = PROJECT_ROOT / "scripts" / "popgen" / "polymorphism_analysis" / "plot_afs_by_class.py"
    output:
        tsv = SFS_DIR / "afs_by_class.tsv",
        pdf = FIGURES_DIR / "snp_popgen" / "afs_by_class.pdf",
        png = FIGURES_DIR / "snp_popgen" / "afs_by_class.png"
    params:
        group1_str = ",".join(GROUP1_SAMPLES),
        outgroup_str = ",".join(GROUP2_SAMPLES)
    log:
        SELECTION_LOG_DIR / "afs_by_class.log"
    conda: "../envs/popgen.yaml"
    shell:
        """
        python {input.script} \
            --snpeff_vcf {input.snpeff_vcf} \
            --genotype_matrix {input.genotype_matrix} \
            --group1 {params.group1_str} \
            --outgroup {params.outgroup_str} \
            --output_tsv {output.tsv} \
            --output_pdf {output.pdf} \
            --output_png {output.png} \
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
    conda: "../envs/popgen.yaml"
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
    conda: "../envs/popgen.yaml"
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
    conda: "../envs/popgen.yaml"
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
        pdf = FIGURES_DIR / "snp_popgen" / "pca.pdf",
        png = FIGURES_DIR / "snp_popgen" / "pca.png"
    params:
        group1_str = ",".join(GROUP1_SAMPLES),
        group2_str = ",".join(GROUP2_SAMPLES)
    log:
        SELECTION_LOG_DIR / "plot_pca.log"
    conda: "../envs/popgen.yaml"
    shell:
        """
        python {PROJECT_ROOT}/scripts/popgen/pca/plot_PCA.py \
            --eigenval {input.eigenval} \
            --eigenvec {input.eigenvec} \
            --output_pdf {output.pdf} \
            --output_png {output.png} \
            --group1 {params.group1_str} \
            --group2 {params.group2_str} \
            2> {log}
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
    conda: "../envs/popgen.yaml"
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
    conda: "../envs/popgen.yaml"
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
        pdf = FIGURES_DIR / "snp_popgen" / "ld_decay.pdf",
        png = FIGURES_DIR / "snp_popgen" / "ld_decay.png"
    log:
        SELECTION_LOG_DIR / "plot_ld.log"
    conda: "../envs/popgen.yaml"
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

rule compare_recombination_introners_all:
    """
    Gene-centric recombination comparison: all introners vs non-introner windows.
    Produces gene_exonic_introners_all_5kb_updated_{summary,windows}.tsv + pdf.
    """
    input:
        gtf = ANNOTATIONS_DIR / "CCMP1545.gtf",
        introner_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv"
    output:
        summary_tsv = RECOMB_DIR / "gene_exonic_introners_all_5kb_updated_summary.tsv",
        windows_tsv = RECOMB_DIR / "gene_exonic_introners_all_5kb_updated_windows.tsv",
        pdf = RECOMB_DIR / "gene_exonic_introners_all_5kb_updated.pdf",
        png = RECOMB_DIR / "gene_exonic_introners_all_5kb_updated.png"
    params:
        pyrho_dir = PYRHO_CCMP1545,
        sample_name = "CCMP1545",
        introner_type = "all",
        window_size = 5000,
        merge_distance = 10000,
        output_prefix = RECOMB_DIR / "gene_exonic_introners_all_5kb_updated"
    log:
        SELECTION_LOG_DIR / "recombination_introners_all.log"
    conda: "../envs/popgen.yaml"
    shell:
        """
        mkdir -p {RECOMB_DIR}

        python {PROJECT_ROOT}/scripts/popgen/recombination_analysis/compare_recombination_gene_exonic_introners.py \
            --pyrho_dir {params.pyrho_dir} \
            --gtf_file {input.gtf} \
            --introner_matrix {input.introner_matrix} \
            --sample_name {params.sample_name} \
            --introner_type {params.introner_type} \
            --window_size {params.window_size} \
            --merge_distance {params.merge_distance} \
            --output_prefix {params.output_prefix} \
            2> {log}
        """


rule compare_recombination_introners_polymorphic:
    """
    Gene-centric recombination comparison: polymorphic introners only.
    Produces gene_exonic_polymorphic_introners_5kb_updated_{summary,windows}.tsv + pdf.
    """
    input:
        gtf = ANNOTATIONS_DIR / "CCMP1545.gtf",
        introner_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv"
    output:
        summary_tsv = RECOMB_DIR / "gene_exonic_polymorphic_introners_5kb_updated_summary.tsv",
        windows_tsv = RECOMB_DIR / "gene_exonic_polymorphic_introners_5kb_updated_windows.tsv",
        pdf = RECOMB_DIR / "gene_exonic_polymorphic_introners_5kb_updated.pdf",
        png = RECOMB_DIR / "gene_exonic_polymorphic_introners_5kb_updated.png"
    params:
        pyrho_dir = PYRHO_CCMP1545,
        sample_name = "CCMP1545",
        introner_type = "polymorphic",
        window_size = 5000,
        merge_distance = 10000,
        output_prefix = RECOMB_DIR / "gene_exonic_polymorphic_introners_5kb_updated"
    log:
        SELECTION_LOG_DIR / "recombination_introners_polymorphic.log"
    conda: "../envs/popgen.yaml"
    shell:
        """
        mkdir -p {RECOMB_DIR}

        python {PROJECT_ROOT}/scripts/popgen/recombination_analysis/compare_recombination_gene_exonic_introners.py \
            --pyrho_dir {params.pyrho_dir} \
            --gtf_file {input.gtf} \
            --introner_matrix {input.introner_matrix} \
            --sample_name {params.sample_name} \
            --introner_type {params.introner_type} \
            --window_size {params.window_size} \
            --merge_distance {params.merge_distance} \
            --output_prefix {params.output_prefix} \
            2> {log}
        """


rule compare_recombination_frequency_based:
    """
    Frequency-based recombination comparison: recent gain vs recent loss vs non-introner.
    Produces gene_exonic_frequency_based_5kb_updated_{summary,windows}.tsv + pdf.
    """
    input:
        gtf = ANNOTATIONS_DIR / "CCMP1545.gtf",
        introner_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv"
    output:
        summary_tsv = RECOMB_DIR / "gene_exonic_frequency_based_5kb_updated_summary.tsv",
        windows_tsv = RECOMB_DIR / "gene_exonic_frequency_based_5kb_updated_windows.tsv",
        pdf = RECOMB_DIR / "gene_exonic_frequency_based_5kb_updated.pdf",
        png = RECOMB_DIR / "gene_exonic_frequency_based_5kb_updated.png"
    params:
        pyrho_dir = PYRHO_CCMP1545,
        exclude_samples = ",".join(GROUP2_SAMPLES),
        window_size = 5000,
        merge_distance = 10000,
        output_prefix = RECOMB_DIR / "gene_exonic_frequency_based_5kb_updated"
    log:
        SELECTION_LOG_DIR / "recombination_frequency_based.log"
    conda: "../envs/popgen.yaml"
    shell:
        """
        mkdir -p {RECOMB_DIR}

        python {PROJECT_ROOT}/scripts/popgen/recombination_analysis/compare_recombination_gene_exonic_frequency_based.py \
            --pyrho_dir {params.pyrho_dir} \
            --gtf_file {input.gtf} \
            --introner_matrix {input.introner_matrix} \
            --exclude_samples {params.exclude_samples} \
            --window_size {params.window_size} \
            --merge_distance {params.merge_distance} \
            --output_prefix {params.output_prefix} \
            2> {log}
        """


rule compare_recombination_group2:
    """
    Recombination at positions where Group 2 (RCC1749/RCC3052) has introners,
    measured in RCC1749's recombination landscape.

    Uses the RCC1749 GTF (contigs named RCC1749#0#intronerless_contig_X) and
    RCC1749 pyrho map together — both are in RCC1749 frame and align with the
    G2-sample contig names in the genotype matrix.
    """
    input:
        gtf = ANNOTATIONS_DIR / "RCC1749.gtf",
        introner_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv"
    output:
        summary_tsv = RECOMB_DIR / "gene_exonic_group2_introners_5kb_updated_summary.tsv",
        windows_tsv = RECOMB_DIR / "gene_exonic_group2_introners_5kb_updated_windows.tsv",
        pdf = RECOMB_DIR / "gene_exonic_group2_introners_5kb_updated.pdf",
        png = RECOMB_DIR / "gene_exonic_group2_introners_5kb_updated.png"
    params:
        pyrho_dir = PYRHO_RCC1749,
        window_size = 5000,
        merge_distance = 10000,
        output_prefix = RECOMB_DIR / "gene_exonic_group2_introners_5kb_updated"
    log:
        SELECTION_LOG_DIR / "recombination_group2.log"
    conda: "../envs/popgen.yaml"
    shell:
        """
        mkdir -p {RECOMB_DIR}

        python {PROJECT_ROOT}/scripts/popgen/recombination_analysis/compare_recombination_gene_exonic_group2_introners.py \
            --pyrho_dir {params.pyrho_dir} \
            --gtf_file {input.gtf} \
            --introner_matrix {input.introner_matrix} \
            --window_size {params.window_size} \
            --merge_distance {params.merge_distance} \
            --output_prefix {params.output_prefix} \
            2> {log}
        """


rule plot_recombination_boxplots:
    """
    Generate boxplots comparing recombination rates by introner category:
    non-introner, all introners, polymorphic, recent gain, recent loss, plus Group 2.
    """
    input:
        all_windows = RECOMB_DIR / "gene_exonic_introners_all_5kb_updated_windows.tsv",
        poly_windows = RECOMB_DIR / "gene_exonic_polymorphic_introners_5kb_updated_windows.tsv",
        freq_windows = RECOMB_DIR / "gene_exonic_frequency_based_5kb_updated_windows.tsv",
        group2_windows = RECOMB_DIR / "gene_exonic_group2_introners_5kb_updated_windows.tsv"
    output:
        pdf = FIGURES_DIR / "snp_popgen" / "recombination_boxplots.pdf",
        png = FIGURES_DIR / "snp_popgen" / "recombination_boxplots.png"
    params:
        base_dir = RECOMB_DIR,
        output_prefix = FIGURES_DIR / "snp_popgen" / "recombination_boxplots"
    log:
        SELECTION_LOG_DIR / "plot_recombination.log"
    conda: "../envs/popgen.yaml"
    shell:
        """
        mkdir -p {FIGURES_DIR}/snp_popgen

        python {PROJECT_ROOT}/scripts/popgen/recombination_analysis/plot_recombination_boxplots.py \
            --base_dir {params.base_dir} \
            --output {params.output_prefix} \
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
        # PCA
        PCA_DIR / "mpusilla.snps.eigenvec",
        FIGURES_DIR / "snp_popgen" / "pca.pdf",
        # LD
        LD_DIR / "ld_summary_by_distance.tsv",
        FIGURES_DIR / "snp_popgen" / "ld_decay.pdf",
        # Recombination
        RECOMB_DIR / "gene_exonic_introners_all_5kb_updated_summary.tsv",
        RECOMB_DIR / "gene_exonic_polymorphic_introners_5kb_updated_summary.tsv",
        RECOMB_DIR / "gene_exonic_frequency_based_5kb_updated_summary.tsv",
        RECOMB_DIR / "gene_exonic_group2_introners_5kb_updated_summary.tsv",
        FIGURES_DIR / "snp_popgen" / "recombination_boxplots.pdf",
        # SFS
        SFS_DIR / "afs_by_class.tsv",
        FIGURES_DIR / "snp_popgen" / "afs_by_class.pdf",


rule basic_popgen_only:
    """
    Target: Basic population genetics statistics only.
    """
    input:
        SELECTION_DIR / "tajimas_d.tsv",
        SFS_DIR / "afs_by_class.tsv",
        FIGURES_DIR / "snp_popgen" / "afs_by_class.pdf"


rule pca_only:
    """
    Target: PCA analysis only.
    """
    input:
        PCA_DIR / "mpusilla.snps.eigenvec",
        FIGURES_DIR / "snp_popgen" / "pca.pdf"


rule ld_only:
    """
    Target: LD analysis only.
    """
    input:
        LD_DIR / "ld_summary_by_distance.tsv",
        FIGURES_DIR / "snp_popgen" / "ld_decay.pdf"


rule recombination_only:
    """
    Target: Recombination analysis only.
    """
    input:
        RECOMB_DIR / "gene_exonic_introners_all_5kb_updated_summary.tsv",
        RECOMB_DIR / "gene_exonic_polymorphic_introners_5kb_updated_summary.tsv",
        RECOMB_DIR / "gene_exonic_frequency_based_5kb_updated_summary.tsv",
        RECOMB_DIR / "gene_exonic_group2_introners_5kb_updated_summary.tsv",
        FIGURES_DIR / "snp_popgen" / "recombination_boxplots.pdf"
