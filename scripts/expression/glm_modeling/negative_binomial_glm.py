#!/usr/bin/env python3
"""
Negative Binomial GLM: Gene Expression Prediction with Library Size Normalization

Reads a pre-merged CSV (from prepare_expression_data inline rule) containing:
  gene_id, replicate, raw_count, library_size, log_library_size,
  strain, introner_gain, introner_loss, baseline_introner_count,
  log_cds_length, n_isoforms

Usage (Snakemake):
    python negative_binomial_glm.py \
        --input CSV --output TXT --coefficients CSV --diagnostics PDF \
        [--reference CCMP1545]
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import statsmodels.api as sm
import statsmodels.formula.api as smf
from patsy.contrasts import Treatment
from scipy import stats
import warnings

warnings.filterwarnings('ignore')


def main():
    parser = argparse.ArgumentParser(
        description="Negative Binomial GLM for gene expression prediction"
    )
    parser.add_argument('--input', required=True, help='Input merged CSV')
    parser.add_argument('--output', required=True, help='Output summary TXT')
    parser.add_argument('--coefficients', required=True, help='Output coefficients CSV')
    parser.add_argument('--diagnostics', required=True, help='Output diagnostic PDF')
    parser.add_argument('--reference', default='CCMP1545',
                        help='Reference strain (default: CCMP1545)')

    args = parser.parse_args()

    # Load pre-merged data
    merged = pd.read_csv(args.input)
    print(f"Loaded {len(merged)} observations")
    print(f"  Unique genes: {merged['gene_id'].nunique()}")
    print(f"  Strains: {merged['strain'].unique().tolist()}")

    # Ensure log_cds_length exists
    if 'log_cds_length' not in merged.columns and 'cds_length' in merged.columns:
        merged['log_cds_length'] = np.log(merged['cds_length'])

    # Build formula - only include available columns
    ref = args.reference
    formula_parts = [f"C(strain, Treatment(reference='{ref}'))"]

    for var in ['introner_gain', 'introner_loss', 'baseline_introner_count',
                'log_cds_length', 'n_isoforms']:
        if var in merged.columns:
            formula_parts.append(var)

    # Include GC_content only if available
    if 'GC_content' in merged.columns:
        formula_parts.append('GC_content')

    formula = 'raw_count ~ ' + ' + '.join(formula_parts)
    print(f"Formula: {formula}")
    print(f"Offset: log_library_size")

    # Fit model
    try:
        model = smf.glm(
            formula,
            data=merged,
            family=sm.families.NegativeBinomial(),
            offset=merged['log_library_size']
        ).fit()
        print(model.summary())
    except Exception as e:
        print(f"Error fitting model: {e}")
        # Write error to output
        with open(args.output, 'w') as f:
            f.write(f"Model fitting failed: {e}\n")
        # Create empty coefficient file
        pd.DataFrame(columns=['variable', 'coefficient', 'std_error',
                               'z_value', 'p_value']).to_csv(args.coefficients, index=False)
        # Create empty diagnostics PDF
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, f'Model fitting failed:\n{e}',
                ha='center', va='center', transform=ax.transAxes)
        fig.savefig(args.diagnostics)
        plt.close()
        return

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
        f.write("Negative Binomial GLM: Gene Expression with Library Size Normalization\n")
        f.write("=" * 80 + "\n\n")

        f.write("MODEL SPECIFICATION\n")
        f.write("-" * 80 + "\n")
        f.write(f"Formula: {formula}\n")
        f.write(f"Family: Negative Binomial\n")
        f.write(f"Offset: log_library_size\n")
        f.write(f"Reference strain: {ref}\n\n")

        f.write("DATASET SUMMARY\n")
        f.write("-" * 80 + "\n")
        f.write(f"Total observations: {len(merged)} (gene x replicate)\n")
        f.write(f"Unique genes: {merged['gene_id'].nunique()}\n")
        f.write(f"Strains: {', '.join(sorted(merged['strain'].unique().tolist()))}\n\n")

        f.write("Raw Count (Response Variable):\n")
        f.write(f"  Mean: {merged['raw_count'].mean():.2f}\n")
        f.write(f"  Median: {merged['raw_count'].median():.2f}\n")
        f.write(f"  Range: {merged['raw_count'].min():.0f} - {merged['raw_count'].max():.0f}\n")
        f.write(f"  Zero counts: {(merged['raw_count'] == 0).sum()} "
                f"({(merged['raw_count'] == 0).mean()*100:.1f}%)\n\n")

        f.write("Library Size:\n")
        f.write(f"  Mean: {merged['library_size'].mean():,.0f}\n")
        f.write(f"  Range: {merged['library_size'].min():,.0f} - "
                f"{merged['library_size'].max():,.0f}\n\n")

        f.write("=" * 80 + "\n")
        f.write("MODEL RESULTS\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"Log-likelihood: {model.llf:.4f}\n")
        f.write(f"AIC: {model.aic:.2f}\n")
        f.write(f"BIC: {model.bic:.2f}\n")
        f.write(f"Deviance: {model.deviance:.2f}\n")
        f.write(f"Dispersion: {model.scale:.4f}\n\n")

        f.write("Coefficients:\n\n")
        f.write(str(model.summary()))
        f.write("\n\n")

        f.write("=" * 80 + "\n")
        f.write("INTERPRETATION\n")
        f.write("=" * 80 + "\n\n")
        f.write("Coefficients represent log-fold changes in expected gene expression.\n")
        f.write("Exponentiate to obtain multiplicative effects on expected counts.\n")
        f.write("Library size offset normalizes for sequencing depth differences.\n")

    print(f"Summary saved to {args.output}")

    # Generate diagnostic plots
    pearson_resid = model.resid_pearson
    deviance_resid = model.resid_deviance
    fitted_values = model.fittedvalues

    with PdfPages(args.diagnostics) as pdf:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        # 1. Residuals vs Fitted
        ax = axes[0, 0]
        ax.scatter(fitted_values, deviance_resid, alpha=0.1, s=5)
        ax.axhline(y=0, color='r', linestyle='--', linewidth=2)
        ax.set_xlabel('Fitted values')
        ax.set_ylabel('Deviance residuals')
        ax.set_title('Residuals vs Fitted')

        # 2. Q-Q plot
        ax = axes[0, 1]
        stats.probplot(deviance_resid, dist="norm", plot=ax)
        ax.set_title('Normal Q-Q (Deviance Residuals)')

        # 3. Scale-Location
        ax = axes[1, 0]
        resid_abs_sqrt = np.sqrt(np.abs(deviance_resid))
        ax.scatter(fitted_values, resid_abs_sqrt, alpha=0.1, s=5)
        ax.set_xlabel('Fitted values')
        ax.set_ylabel('sqrt(|Deviance residuals|)')
        ax.set_title('Scale-Location')

        # 4. Histogram
        ax = axes[1, 1]
        ax.hist(deviance_resid, bins=50, color='skyblue', edgecolor='black', alpha=0.7)
        ax.axvline(x=0, color='r', linestyle='--', linewidth=2)
        ax.set_xlabel('Deviance residuals')
        ax.set_ylabel('Frequency')
        ax.set_title('Histogram of Residuals')

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

    print(f"Diagnostic plots saved to {args.diagnostics}")


if __name__ == '__main__':
    main()
