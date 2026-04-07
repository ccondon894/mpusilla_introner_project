#!/usr/bin/env python3
"""
Poisson GLM: Gene Length and Introner Interactions

Tests whether CDS length moderates the effect of introners on isoform diversity.

Usage (Snakemake):
    python poisson_glm_gene_length.py \
        --input CSV --output TXT --coefficients CSV [--reference CCMP1545]
"""

import argparse
import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
from patsy.contrasts import Treatment
import warnings

warnings.filterwarnings('ignore')


def main():
    parser = argparse.ArgumentParser(
        description="Poisson GLM with gene length interaction effects"
    )
    parser.add_argument('--input', required=True, help='Input CSV data')
    parser.add_argument('--output', required=True, help='Output summary TXT')
    parser.add_argument('--coefficients', required=True, help='Output coefficients CSV')
    parser.add_argument('--reference', default='CCMP1545',
                        help='Reference strain (default: CCMP1545)')

    args = parser.parse_args()

    # Load data
    data = pd.read_csv(args.input)
    print(f"Dataset: {len(data)} observations, {data['gene_id'].nunique()} unique genes")
    print(f"Strains: {data['strain'].unique().tolist()}")

    # Filter extreme outliers
    n_before = len(data)
    data = data[data['n_isoforms'] <= 20]
    print(f"Removed {n_before - len(data)} extreme outliers (>20 isoforms)")

    # Create log-transformed variables if not present
    if 'log_expression' not in data.columns:
        if 'mean_expression' in data.columns:
            data['log_expression'] = np.log(data['mean_expression'] + 1)

    if 'log_cds_length' not in data.columns:
        data['log_cds_length'] = np.log(data['cds_length'])

    # Build formula
    ref = args.reference
    formula_parts = [
        f"C(strain, Treatment(reference='{ref}'))",
        'introner_loss', 'introner_gain', 'baseline_introner_count',
        'log_cds_length'
    ]

    if 'log_expression' in data.columns:
        formula_parts.append('log_expression')

    formula = 'n_isoforms ~ ' + ' + '.join(formula_parts)
    print(f"Formula: {formula}")

    # Fit model
    model = smf.glm(
        formula,
        data=data,
        family=sm.families.Poisson()
    ).fit()

    print(model.summary())

    # Save coefficients
    coef_table = pd.DataFrame({
        'variable': model.params.index,
        'coefficient': model.params.values,
        'std_error': model.bse.values,
        'z_value': model.tvalues.values,
        'p_value': model.pvalues.values,
        'ci_lower': model.conf_int()[0].values,
        'ci_upper': model.conf_int()[1].values
    })
    coef_table.to_csv(args.coefficients, index=False)
    print(f"Coefficients saved to {args.coefficients}")

    # Save summary report
    with open(args.output, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("Poisson GLM: CDS Length and Introner Interactions\n")
        f.write("=" * 80 + "\n\n")

        f.write("MODEL SPECIFICATION\n")
        f.write("-" * 80 + "\n")
        f.write(f"Formula: {formula}\n")
        f.write("Family: Poisson (log link)\n")
        f.write(f"Reference strain: {ref}\n\n")

        f.write("DATASET SUMMARY\n")
        f.write("-" * 80 + "\n")
        f.write(f"Total observations: {len(data)}\n")
        f.write(f"Unique genes: {data['gene_id'].nunique()}\n")
        f.write(f"Strains: {', '.join(data['strain'].unique().tolist())}\n\n")

        f.write(f"n_isoforms: mean={data['n_isoforms'].mean():.3f}, "
                f"median={data['n_isoforms'].median():.0f}\n")
        f.write(f"CDS Length: mean={data['cds_length'].mean():.0f} bp, "
                f"median={data['cds_length'].median():.0f} bp\n")
        f.write(f"Baseline introner count: mean={data['baseline_introner_count'].mean():.3f}\n\n")

        f.write("=" * 80 + "\n")
        f.write("MODEL RESULTS\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"AIC: {model.aic:.2f}\n")
        f.write(f"BIC: {model.bic:.2f}\n")
        f.write(f"Log-likelihood: {model.llf:.4f}\n")
        f.write(f"Deviance: {model.deviance:.2f}\n\n")

        f.write("Coefficients:\n\n")
        f.write(str(model.summary()))
        f.write("\n")

    print(f"Summary saved to {args.output}")


if __name__ == '__main__':
    main()
