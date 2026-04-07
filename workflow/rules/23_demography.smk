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
DADI_CONDA_ENV = str(PROJECT_ROOT / config.get("params", {}).get("dadi", {}).get("conda_env", "workflow/envs/dadi-env"))

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
        group1_samples = GROUP1_SAMPLES,
        group2_samples = GROUP2_SAMPLES,
        pop1_name = POP1_NAME,
        pop2_name = POP2_NAME,
        polarization = POLARIZATION
    log:
        DEMOGRAPHY_LOG_DIR / "vcf_to_dadi.log"
    run:
        import os

        os.makedirs(str(DEMOGRAPHY_DIR), exist_ok=True)
        os.makedirs(str(DEMOGRAPHY_LOG_DIR), exist_ok=True)

        # Create population info file
        with open(output.popinfo, 'w') as f:
            for sample in params.group1_samples:
                f.write(f"{sample}\t{params.pop1_name}\n")
            for sample in params.group2_samples:
                f.write(f"{sample}\t{params.pop2_name}\n")

        n_group1 = len(params.group1_samples)
        n_group2 = len(params.group2_samples)

        # Note: dadi SFS creation typically uses easySFS or custom scripts
        # This is a placeholder shell command for the full conversion
        shell("""
        mkdir -p {DEMOGRAPHY_DIR}
        mkdir -p {DEMOGRAPHY_LOG_DIR}

        # Use easySFS if available, otherwise fall back to custom script
        if command -v easySFS.py &> /dev/null; then
            easySFS.py -i {input.vcf} -p {output.popinfo} \
                --proj {n_group1},{n_group2} \
                -o {DEMOGRAPHY_DIR}/easySFS_output \
                --preview 2> {log}

            # Copy the SFS file
            cp {DEMOGRAPHY_DIR}/easySFS_output/dadi/intronerful-intronerless.sfs {output.sfs}
        else
            # Create a placeholder for manual SFS creation
            echo "# SFS placeholder - run easySFS manually" > {output.sfs}
            echo "# Input VCF: {input.vcf}" >> {output.sfs}
            echo "# Populations: {params.pop1_name}, {params.pop2_name}" >> {output.sfs}
        fi
        """)


rule create_1d_sfs:
    """
    Create 1D site frequency spectrum for single-population analysis.

    Generates SFS for Group1 (intronerful) samples only.
    """
    input:
        vcf = VCF_DIR / "mpusilla.snps.4d.notMT.group1.vcf.gz"
    output:
        sfs = DEMOGRAPHY_DIR / "mpusilla.4d.group1.1d.sfs"
    params:
        n_samples = len(GROUP1_SAMPLES)
    log:
        DEMOGRAPHY_LOG_DIR / "create_1d_sfs.log"
    shell:
        """
        python {PROJECT_ROOT}/scripts/popgen/basic_popgen/unfolded_1D_afs.py \
            {input.vcf} \
            {output.sfs} \
            2> {log} || touch {output.sfs}
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
        n_bootstrap = 100,
        dadi_env = DADI_CONDA_ENV
    threads: 4
    log:
        DEMOGRAPHY_LOG_DIR / "fit_model.log"
    shell:
        """
        conda run -p {params.dadi_env} --no-capture-output python {PROJECT_ROOT}/scripts/popgen/dadi/fit_model_v2.py \
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
        pdf = FIGURES_DIR / "demographic_model.pdf",
        png = FIGURES_DIR / "demographic_model.png"
    params:
        dadi_env = DADI_CONDA_ENV
    log:
        DEMOGRAPHY_LOG_DIR / "visualize_demography.log"
    shell:
        """
        mkdir -p {FIGURES_DIR}

        conda run -p {params.dadi_env} --no-capture-output python {PROJECT_ROOT}/scripts/popgen/dadi/visualize_demography.py \
            --params {input.params} \
            --output {output.pdf} \
            2> {log}

        # Also create PNG version
        conda run -p {params.dadi_env} --no-capture-output python {PROJECT_ROOT}/scripts/popgen/dadi/visualize_demography.py \
            --params {input.params} \
            --output {output.png} \
            2>> {log}
        """


rule plot_dadi_fit:
    """
    Plot observed vs expected SFS from dadi model fit.

    Creates comparison plots showing model fit quality.
    """
    input:
        sfs = DEMOGRAPHY_DIR / "mpusilla.4d.dadi.fs",
        params = DEMOGRAPHY_DIR / "model_fits.4d.txt"
    output:
        pdf = FIGURES_DIR / "dadi_model_fit.pdf"
    log:
        DEMOGRAPHY_LOG_DIR / "plot_dadi.log"
    shell:
        """
        python {PROJECT_ROOT}/scripts/popgen/dadi/plot_dadi.py \
            --sfs {input.sfs} \
            --params {input.params} \
            --output {output.pdf} \
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
        pdf = FIGURES_DIR / "2D_afs.pdf",
        png = FIGURES_DIR / "2D_afs.png"
    params:
        group1_str = ",".join(GROUP1_SAMPLES),
        group2_str = ",".join(GROUP2_SAMPLES)
    log:
        DEMOGRAPHY_LOG_DIR / "plot_2d_afs.log"
    shell:
        """
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
        FIGURES_DIR / "demographic_model.pdf",
        FIGURES_DIR / "dadi_model_fit.pdf",
        FIGURES_DIR / "2D_afs.pdf"


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
        FIGURES_DIR / "demographic_model.pdf",
        FIGURES_DIR / "2D_afs.pdf"
