# 23_demography.smk - Demographic Modeling with moments

import os
from pathlib import Path


# Output directories
SNP_DIR = RESULTS / "snp_popgen"
VCF_DIR = SNP_DIR / "vcf"
DEMOGRAPHY_DIR = SNP_DIR / "demography"
DEMOGRAPHY_LOG_DIR = SNP_DIR / "logs" / "demography"

# Population definitions for moments
# Group1 = intronerful, Group2 = intronerless
POP1_NAME = "intronerful"
POP2_NAME = "intronerless"


# Demographic model fitting

rule fit_demographic_model:
    """Fit split-migration demographic model using moments."""
    input:
        vcf = VCF_DIR / "mpusilla.snps.4d.notMT.vcf.gz"
    output:
        popinfo = DEMOGRAPHY_DIR / "popinfo.txt",
        sfs = DEMOGRAPHY_DIR / "observed.2d.fs",
        params = DEMOGRAPHY_DIR / "model_fits.4d.txt",
        bootstrap = DEMOGRAPHY_DIR / "model_fits.4d.bootstrap.txt",
        summary = DEMOGRAPHY_DIR / "bootstrap_summary.tsv"
    params:
        group1_samples = " ".join(GROUP1_SAMPLES),
        group2_samples = " ".join(GROUP2_SAMPLES),
        n_group1 = len(GROUP1_SAMPLES),
        n_group2 = len(GROUP2_SAMPLES),
        pop1_name = POP1_NAME,
        pop2_name = POP2_NAME,
        n_optimizations = 20,
        n_bootstrap = 100,
        chunk_size = 250000,
        maxiter = 10000
    threads: 4
    log:
        DEMOGRAPHY_LOG_DIR / "fit_moments.log"
    conda: "../envs/moments.yaml"
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

        python {PROJECT_ROOT}/scripts/popgen/moments/fit_model.py \
            --vcf {input.vcf} \
            --popinfo {output.popinfo} \
            --outdir {DEMOGRAPHY_DIR} \
            --output {output.params} \
            --bootstrap-output {output.bootstrap} \
            --pop-ids {params.pop1_name} {params.pop2_name} \
            --projections {params.n_group1} {params.n_group2} \
            --n-opt {params.n_optimizations} \
            --n-boot {params.n_bootstrap} \
            --chunk-size {params.chunk_size} \
            --maxiter {params.maxiter} \
            --threads {threads} \
            2> {log}
        """


# Visualization

rule visualize_demography:
    """Create demographic model visualization."""
    input:
        params = DEMOGRAPHY_DIR / "model_fits.4d.txt",
        target_bed = SNP_DIR / "degenotate" / "degeneracy-all-sites.4d.bed.gz"
    output:
        pdf = FIGURES_DIR / "snp_popgen" / "demographic_model.pdf",
        png = FIGURES_DIR / "snp_popgen" / "demographic_model.png"
    log:
        DEMOGRAPHY_LOG_DIR / "visualize_demography.log"
    conda: "../envs/moments.yaml"
    shell:
        """
        mkdir -p {FIGURES_DIR}/snp_popgen

        python {PROJECT_ROOT}/scripts/popgen/dadi/visualize_demography.py \
            --params {input.params} \
            --target-bed {input.target_bed} \
            --mutation-rate 9.8e-10 \
            --ancestral-time-factor 40 \
            --figure-width 6.5 \
            --figure-height 2.35 \
            --font-size 10 \
            --color-guide {PROJECT_ROOT}/master_figure_color_guide.tsv \
            --output_pdf {output.pdf} \
            --output_png {output.png} \
            2> {log}
        """


rule plot_moments_fit:
    """Plot scatter matrix of moments bootstrap parameter estimates with 95% CIs."""
    input:
        bootstrap = DEMOGRAPHY_DIR / "model_fits.4d.bootstrap.txt"
    output:
        pdf = FIGURES_DIR / "snp_popgen" / "moments_model_fit.pdf",
        png = FIGURES_DIR / "snp_popgen" / "moments_model_fit.png"
    log:
        DEMOGRAPHY_LOG_DIR / "plot_moments.log"
    conda: "../envs/moments.yaml"
    shell:
        """
        mkdir -p {FIGURES_DIR}/snp_popgen

        python {PROJECT_ROOT}/scripts/popgen/dadi/plot_dadi.py \
            {input.bootstrap} \
            --output_file {output.pdf} \
            2> {log}
        """


rule plot_2d_afs:
    """Plot 2D allele frequency spectrum."""
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
            --figure-width 6.5 \
            --figure-height 2.7 \
            --font-size 10 \
            --output_pdf {output.pdf} \
            --output_png {output.png} \
            2> {log}
        """


rule demography_complete:
    """Target: Complete demographic analysis."""
    input:
        DEMOGRAPHY_DIR / "observed.2d.fs",
        DEMOGRAPHY_DIR / "model_fits.4d.txt",
        FIGURES_DIR / "snp_popgen" / "demographic_model.pdf",
        FIGURES_DIR / "snp_popgen" / "moments_model_fit.pdf",
        FIGURES_DIR / "snp_popgen" / "2D_afs.pdf"


rule moments_fitting_only:
    """Target: Run moments model fitting only."""
    input:
        DEMOGRAPHY_DIR / "model_fits.4d.txt"


rule demography_plots_only:
    """Target: Generate demographic plots only."""
    input:
        FIGURES_DIR / "snp_popgen" / "demographic_model.pdf",
        FIGURES_DIR / "snp_popgen" / "2D_afs.pdf"
