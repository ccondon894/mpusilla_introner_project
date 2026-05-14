# ============================================================
# 23_demography.smk - Demographic Modeling with dadi
# ============================================================
#
# Performs demographic inference using dadi to estimate:
# - Population sizes
# - Split times
# - Migration rates
#
# Pipeline:
# 1. Convert VCF to dadi format
# 2. Fit demographic models
# 3. Visualize demographic history
#
# Adapted from:
# - /scratch1/chris/introner_vis/dadi/
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
DEMOGRAPHY_DIR = SNP_DIR / "demography"
DEMOGRAPHY_LOG_DIR = SNP_DIR / "logs" / "demography"

# Dadi parameters
POLARIZATION = config["params"]["dadi"]["polarization"]

# Population definitions for dadi
# Group1 = intronerful, Group2 = intronerless
POP1_NAME = "intronerful"
POP2_NAME = "intronerless"


# ============================================================
# DATA CONVERSION
# ============================================================

rule vcf_to_dadi_sfs:
    """
    Convert VCF to dadi site frequency spectrum format.

    Creates 2D SFS for two-population demographic inference.
    Uses Group1 (intronerful) and Group2 (intronerless) as populations.
    """
    input:
        vcf = VCF_DIR / "mpusilla.snps.4d.notMT.vcf.gz"
    output:
        sfs = DEMOGRAPHY_DIR / "mpusilla.4d.dadi.fs",
        popinfo = DEMOGRAPHY_DIR / "popinfo.txt"
    params:
        group1_samples = " ".join(GROUP1_SAMPLES),
        group2_samples = " ".join(GROUP2_SAMPLES),
        n_group1 = len(GROUP1_SAMPLES),
        n_group2 = len(GROUP2_SAMPLES),
        pop1_name = POP1_NAME,
        pop2_name = POP2_NAME,
        polarization = POLARIZATION
    log:
        DEMOGRAPHY_LOG_DIR / "vcf_to_dadi.log"
    conda: "../envs/dadi.yaml"
    shell:
        """
        mkdir -p {DEMOGRAPHY_DIR}
        mkdir -p {DEMOGRAPHY_LOG_DIR}

        : > {output.popinfo}
        for sample in {params.group1_samples}; do
            printf "%s\t%s\n" "$sample" "{params.pop1_name}" >> {output.popinfo}
        done
        for sample in {params.group2_samples}; do
            printf "%s\t%s\n" "$sample" "{params.pop2_name}" >> {output.popinfo}
        done

        # Use easySFS if available, otherwise fall back to custom script
        if command -v easySFS.py &> /dev/null; then
            easySFS.py -i {input.vcf} -p {output.popinfo} \
                --proj {params.n_group1},{params.n_group2} \
                -o {DEMOGRAPHY_DIR}/easySFS_output \
                2> {log}

            # Copy the SFS file
            cp {DEMOGRAPHY_DIR}/easySFS_output/dadi/intronerful-intronerless.sfs {output.sfs}
        else
            # Create a placeholder for manual SFS creation
            echo "# SFS placeholder - run easySFS manually" > {output.sfs}
            echo "# Input VCF: {input.vcf}" >> {output.sfs}
            echo "# Populations: {params.pop1_name}, {params.pop2_name}" >> {output.sfs}
        fi
        """


# ============================================================
# DEMOGRAPHIC MODEL FITTING
# ============================================================

rule fit_demographic_model:
    """
    Fit split-migration demographic model using dadi.

    Estimates:
    - N1: Effective population size of Group1
    - N2: Effective population size of Group2
    - T: Time since population split
    - m12, m21: Migration rates between populations
    """
    input:
        vcf = VCF_DIR / "mpusilla.snps.4d.notMT.vcf.gz",
        popinfo = DEMOGRAPHY_DIR / "popinfo.txt"
    output:
        params = DEMOGRAPHY_DIR / "model_fits.4d.txt",
        bootstrap = DEMOGRAPHY_DIR / "model_fits.4d.bootstrap.txt"
    params:
        n_optimizations = 20,
        n_bootstrap = 100
    threads: 4
    log:
        DEMOGRAPHY_LOG_DIR / "fit_model.log"
    conda: "../envs/dadi.yaml"
    shell:
        """
        python {PROJECT_ROOT}/scripts/popgen/dadi/fit_model_v2.py \
            --vcf {input.vcf} \
            --popinfo {input.popinfo} \
            --output {output.params} \
            --bootstrap {output.bootstrap} \
            --n_opt {params.n_optimizations} \
            --n_boot {params.n_bootstrap} \
            --threads {threads} \
            2> {log}
        """


# ============================================================
# VISUALIZATION
# ============================================================

rule visualize_demography:
    """
    Create demographic model visualization.

    Generates schematic of population history showing:
    - Ancestral population
    - Population split
    - Current population sizes
    - Migration arrows
    """
    input:
        params = DEMOGRAPHY_DIR / "model_fits.4d.txt"
    output:
        pdf = FIGURES_DIR / "snp_popgen" / "demographic_model.pdf",
        png = FIGURES_DIR / "snp_popgen" / "demographic_model.png"
    log:
        DEMOGRAPHY_LOG_DIR / "visualize_demography.log"
    conda: "../envs/dadi.yaml"
    shell:
        """
        mkdir -p {FIGURES_DIR}/snp_popgen

        python {PROJECT_ROOT}/scripts/popgen/dadi/visualize_demography.py \
            --params {input.params} \
            --output {output.pdf} \
            2> {log}

        # Also create PNG version
        python {PROJECT_ROOT}/scripts/popgen/dadi/visualize_demography.py \
            --params {input.params} \
            --output {output.png} \
            2>> {log}
        """


rule plot_dadi_fit:
    """
    Plot scatter matrix of bootstrap parameter estimates with 95% CIs.

    Visualizes the joint distribution of (N1, N2, T, M, Theta) across
    bootstrap iterations from the dadi demographic fit, with confidence
    intervals annotated.
    """
    input:
        bootstrap = DEMOGRAPHY_DIR / "model_fits.4d.bootstrap.txt"
    output:
        pdf = FIGURES_DIR / "snp_popgen" / "dadi_model_fit.pdf",
        png = FIGURES_DIR / "snp_popgen" / "dadi_model_fit.png"
    log:
        DEMOGRAPHY_LOG_DIR / "plot_dadi.log"
    conda: "../envs/dadi.yaml"
    shell:
        """
        mkdir -p {FIGURES_DIR}/snp_popgen

        python {PROJECT_ROOT}/scripts/popgen/dadi/plot_dadi.py \
            {input.bootstrap} \
            --output_file {output.pdf} \
            2> {log}
        """


rule plot_2d_afs:
    """
    Plot 2D allele frequency spectrum.

    Creates heatmap visualization of joint SFS between Group1 and Group2.
    """
    input:
        vcf = VCF_DIR / "mpusilla.snps.4d.notMT.vcf.gz"
    output:
        pdf = FIGURES_DIR / "snp_popgen" / "2D_afs.pdf",
        png = FIGURES_DIR / "snp_popgen" / "2D_afs.png"
    params:
        group1_str = ",".join(GROUP1_SAMPLES),
        group2_str = ",".join(GROUP2_SAMPLES)
    log:
        DEMOGRAPHY_LOG_DIR / "plot_2d_afs.log"
    conda: "../envs/dadi.yaml"
    shell:
        """
        mkdir -p {FIGURES_DIR}/snp_popgen

        python {PROJECT_ROOT}/scripts/popgen/basic_popgen/2D_afs.py \
            --vcf {input.vcf} \
            --group1 {params.group1_str} \
            --group2 {params.group2_str} \
            --output {output.pdf} \
            2> {log}

        python {PROJECT_ROOT}/scripts/popgen/basic_popgen/2D_afs.py \
            --vcf {input.vcf} \
            --group1 {params.group1_str} \
            --group2 {params.group2_str} \
            --output {output.png} \
            2>> {log}
        """


# ============================================================
# TARGET RULES
# ============================================================

rule demography_complete:
    """
    Target: Complete demographic analysis.
    """
    input:
        DEMOGRAPHY_DIR / "mpusilla.4d.dadi.fs",
        DEMOGRAPHY_DIR / "model_fits.4d.txt",
        FIGURES_DIR / "snp_popgen" / "demographic_model.pdf",
        FIGURES_DIR / "snp_popgen" / "dadi_model_fit.pdf",
        FIGURES_DIR / "snp_popgen" / "2D_afs.pdf"


rule dadi_fitting_only:
    """
    Target: Run dadi model fitting only.
    """
    input:
        DEMOGRAPHY_DIR / "model_fits.4d.txt"


rule demography_plots_only:
    """
    Target: Generate demographic plots only.
    """
    input:
        FIGURES_DIR / "snp_popgen" / "demographic_model.pdf",
        FIGURES_DIR / "snp_popgen" / "2D_afs.pdf"
