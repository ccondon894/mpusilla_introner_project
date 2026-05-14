# ============================================================
# 32_glm_modeling.smk - Statistical GLM Modeling
# ============================================================
#
# Performs statistical modeling to test introner effects:
# 1. Poisson GLM for isoform diversity (count response)
# 2. Negative binomial GLM for expression levels
# 3. Gene length interaction effects
#
# Research Questions:
# - Do introner gain/loss predict number of transcript isoforms?
# - How do introners affect gene expression levels?
# - Is there an interaction with gene architecture (CDS length)?
#
# Adapted from:
# - /scratch1/chris/introner-expression-analysis/isoform_analysis/poisson_modeling/
#
# ============================================================

import os
from pathlib import Path

# ============================================================
# CONFIGURATION
# ============================================================

# Output directories
EXPRESSION_DIR = RESULTS / "expression"
GLM_DIR = EXPRESSION_DIR / "glm_modeling"
EXPRESSION_LOG_DIR = EXPRESSION_DIR / "logs"

# Input from isoform analysis
ISOFORM_DIR = EXPRESSION_DIR / "isoform_analysis"
COUNTS_DIR = EXPRESSION_DIR / "counts"


# ============================================================
# POISSON GLM - ISOFORM DIVERSITY
# ============================================================

rule poisson_glm_isoform_diversity:
    """
    Fit Poisson GLM for isoform count (diversity) prediction.

    Model: n_isoforms ~ strain + introner_gain + introner_loss +
                        log_expression + log_cds_length

    Tests whether introner gain/loss predicts number of transcript
    isoforms per gene.
    """
    input:
        data = ISOFORM_DIR / "isoform_introner_data_filtered.csv"
    output:
        results = GLM_DIR / "poisson_glm_summary.txt",
        coefficients = GLM_DIR / "poisson_glm_coefficients.csv",
        diagnostics = GLM_DIR / "poisson_glm_diagnostics.pdf",
        diagnostics_png = GLM_DIR / "poisson_glm_diagnostics.png"
    log:
        EXPRESSION_LOG_DIR / "poisson_glm.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        mkdir -p {GLM_DIR}

        python {PROJECT_ROOT}/scripts/expression/glm_modeling/poisson_GLM_regression.py \
            --input {input.data} \
            --output {output.results} \
            --coefficients {output.coefficients} \
            --diagnostics {output.diagnostics} \
            2> {log}
        """


rule poisson_glm_gene_length:
    """
    Fit Poisson GLM with gene length interaction effects.

    Tests whether CDS length moderates the effect of introners
    on isoform diversity.
    """
    input:
        data = ISOFORM_DIR / "isoform_introner_data_filtered.csv"
    output:
        results = GLM_DIR / "poisson_glm_length_interaction.txt",
        coefficients = GLM_DIR / "poisson_glm_length_coefficients.csv"
    log:
        EXPRESSION_LOG_DIR / "poisson_glm_length.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        python {PROJECT_ROOT}/scripts/expression/glm_modeling/poisson_glm_gene_length.py \
            --input {input.data} \
            --output {output.results} \
            --coefficients {output.coefficients} \
            2> {log}
        """


# ============================================================
# NEGATIVE BINOMIAL GLM - EXPRESSION LEVELS
# ============================================================

rule prepare_expression_data:
    """
    Prepare expression data for negative binomial GLM.

    Merges:
    - Raw counts from featureCounts
    - Library size normalization factors
    - Introner features
    - Gene architecture features (CDS length, GC content)
    """
    input:
        counts = COUNTS_DIR / "merged_counts_matrix.csv",
        isoform_data = ISOFORM_DIR / "isoform_introner_data_filtered.csv"
    output:
        data = GLM_DIR / "expression_glm_data.csv"
    log:
        EXPRESSION_LOG_DIR / "prepare_expression_data.log"
    run:
        import pandas as pd
        import numpy as np

        # Load count matrix
        counts = pd.read_csv(input.counts, index_col=0)

        # Calculate library sizes
        sample_cols = [c for c in counts.columns if c not in ['Length', 'Chr', 'Start', 'End', 'Strand']]
        library_sizes = counts[sample_cols].sum()

        # Melt to long format
        count_long = counts[sample_cols].reset_index().melt(
            id_vars=['Geneid'] if 'Geneid' in counts.index.name else [counts.index.name],
            var_name='replicate',
            value_name='raw_count'
        )
        count_long.columns = ['gene_id', 'replicate', 'raw_count']

        # Add library size
        count_long['library_size'] = count_long['replicate'].map(library_sizes)
        count_long['log_library_size'] = np.log(count_long['library_size'])

        # Map replicate to strain
        strain_map = {}
        for rep in count_long['replicate'].unique():
            if rep.startswith('834'):
                strain_map[rep] = 'CCMP1545'
            elif rep.startswith('1614'):
                strain_map[rep] = 'RCC1614'
            elif rep.startswith('1749'):
                strain_map[rep] = 'RCC1749'
        count_long['strain'] = count_long['replicate'].map(strain_map)

        # Load isoform data for introner features
        isoform_data = pd.read_csv(input.isoform_data)

        # Merge on gene_id and strain
        merged = count_long.merge(
            isoform_data[['gene_id', 'strain', 'introner_gain', 'introner_loss',
                          'baseline_introner_count', 'log_cds_length', 'n_isoforms']].drop_duplicates(),
            on=['gene_id', 'strain'],
            how='inner'
        )

        merged.to_csv(output.data, index=False)


rule negative_binomial_glm:
    """
    Fit Negative Binomial GLM for gene expression.

    Model: raw_count ~ offset(log_library_size) + strain + introner_gain +
                       introner_loss + baseline_introner_count + log_cds_length +
                       n_isoforms + gc_content

    Tests how introners and gene architecture affect expression levels.
    Uses negative binomial to handle overdispersion in count data.
    """
    input:
        data = GLM_DIR / "expression_glm_data.csv"
    output:
        results = GLM_DIR / "negative_binomial_glm_results.txt",
        coefficients = GLM_DIR / "negative_binomial_glm_coefficients.csv",
        diagnostics = GLM_DIR / "negative_binomial_glm_diagnostics.pdf",
        diagnostics_png = GLM_DIR / "negative_binomial_glm_diagnostics.png"
    log:
        EXPRESSION_LOG_DIR / "negative_binomial_glm.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        python {PROJECT_ROOT}/scripts/expression/glm_modeling/negative_binomial_glm.py \
            --input {input.data} \
            --output {output.results} \
            --coefficients {output.coefficients} \
            --diagnostics {output.diagnostics} \
            2> {log}
        """


# ============================================================
# MODEL COMPARISON AND VISUALIZATION
# ============================================================

rule compare_glm_models:
    """
    Compare Poisson vs Negative Binomial model fits.

    Uses AIC/BIC for model selection and likelihood ratio tests.
    """
    input:
        poisson = GLM_DIR / "poisson_glm_summary.txt",
        negbin = GLM_DIR / "negative_binomial_glm_results.txt"
    output:
        comparison = GLM_DIR / "model_comparison.txt",
        plot = FIGURES_DIR / "glm_model_comparison.pdf"
    log:
        EXPRESSION_LOG_DIR / "compare_models.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        python {PROJECT_ROOT}/scripts/expression/glm_modeling/compare_glm_models.py \
            --poisson {input.poisson} \
            --negbin {input.negbin} \
            --output {output.comparison} \
            --plot {output.plot} \
            2> {log}
        """


rule visualize_glm_effects:
    """
    Visualize significant GLM effects.

    Creates coefficient plots and effect size visualizations.
    """
    input:
        poisson_coef = GLM_DIR / "poisson_glm_coefficients.csv",
        negbin_coef = GLM_DIR / "negative_binomial_glm_coefficients.csv"
    output:
        effect_plot = FIGURES_DIR / "glm_effect_sizes.pdf",
        effect_png = FIGURES_DIR / "glm_effect_sizes.png",
        forest_plot = FIGURES_DIR / "glm_forest_plot.pdf",
        forest_png = FIGURES_DIR / "glm_forest_plot.png"
    log:
        EXPRESSION_LOG_DIR / "visualize_glm.log"
    conda: "../envs/glm_modeling.yaml"
    shell:
        """
        python {PROJECT_ROOT}/scripts/expression/glm_modeling/visualize_glm_effects.py \
            --poisson {input.poisson_coef} \
            --negbin {input.negbin_coef} \
            --effect-pdf {output.effect_plot} \
            --effect-png {output.effect_png} \
            --forest-pdf {output.forest_plot} \
            --forest-png {output.forest_png} \
            2> {log}
        """


# ============================================================
# TARGET RULES
# ============================================================

rule glm_modeling_complete:
    """
    Target: Complete GLM statistical modeling.
    """
    input:
        GLM_DIR / "poisson_glm_summary.txt",
        GLM_DIR / "negative_binomial_glm_results.txt",
        GLM_DIR / "model_comparison.txt",
        FIGURES_DIR / "glm_effect_sizes.pdf"


rule poisson_glm_only:
    """
    Target: Poisson GLM analysis only.
    """
    input:
        GLM_DIR / "poisson_glm_summary.txt",
        GLM_DIR / "poisson_glm_length_interaction.txt"


rule expression_glm_only:
    """
    Target: Negative binomial expression GLM only.
    """
    input:
        GLM_DIR / "negative_binomial_glm_results.txt"
