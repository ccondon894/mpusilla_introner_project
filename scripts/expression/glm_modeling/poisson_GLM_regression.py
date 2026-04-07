#!/usr/bin/env python3
"""
Poisson GLM Regression Analysis: Effect of Introners on Isoform Diversity

Model: n_isoforms ~ C(strain) + introner_gain + introner_loss + log_expression
Family: Poisson (log link)

Usage (Snakemake):
    python poisson_GLM_regression.py \
        --input CSV --output TXT --coefficients CSV --diagnostics PDF
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
from statsmodels.genmod.families import Poisson, NegativeBinomial
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy import stats
import warnings

warnings.filterwarnings('ignore')


def main():
    parser = argparse.ArgumentParser(
        description="Poisson GLM for isoform diversity prediction"
    )
    parser.add_argument('--input', required=True, help='Input CSV data')
    parser.add_argument('--output', required=True, help='Output summary TXT')
    parser.add_argument('--coefficients', required=True, help='Output coefficients CSV')
    parser.add_argument('--diagnostics', required=True, help='Output diagnostic plots PDF')

    args = parser.parse_args()

    # Load data
    data = pd.read_csv(args.input)
    print(f"Loaded {len(data)} observations, {data['gene_id'].nunique()} genes")

    # Filter extreme outliers
    n_before = len(data)
    data = data[data['n_isoforms'] <= 20]
    print(f"Removed {n_before - len(data)} extreme outliers (>20 isoforms)")

    # Log-transform expression if available
    if 'log_expression' not in data.columns:
        if 'mean_expression' in data.columns:
            data['log_expression'] = np.log(data['mean_expression'] + 1)

    # Build formula
    has_expr = 'log_expression' in data.columns
    if has_expr:
        formula = 'n_isoforms ~ C(strain) + introner_gain + introner_loss + log_expression'
    else:
        formula = 'n_isoforms ~ C(strain) + introner_gain + introner_loss'

    # Fit Poisson GLM
    print(f"Fitting: {formula}")
    model = smf.glm(formula, data=data, family=Poisson()).fit()
    print(model.summary())

    # Variance-to-mean ratio
    n_isoforms = data['n_isoforms']
    mean_val = n_isoforms.mean()
    var_val = n_isoforms.var()
    vmr = var_val / mean_val

    # Overdispersion check
    dispersion_pearson = model.pearson_chi2 / model.df_resid
    dispersion_deviance = model.deviance / model.df_resid
    overdispersed = dispersion_pearson > 1.5 or dispersion_deviance > 1.5

    # Try Negative Binomial if overdispersed
    best_model = model
    model_name = "Poisson"

    if overdispersed:
        print("Overdispersion detected, fitting Negative Binomial...")
        try:
            nb_model = smf.glm(formula, data=data, family=NegativeBinomial()).fit()
            if nb_model.aic < model.aic - 2:
                best_model = nb_model
                model_name = "Negative Binomial"
                print(f"Negative Binomial preferred (AIC: {nb_model.aic:.1f} vs {model.aic:.1f})")
        except Exception as e:
            print(f"NB fitting failed: {e}, using Poisson")

    # Extract coefficients
    coef_table = pd.DataFrame({
        'variable': best_model.params.index,
        'coefficient': best_model.params.values,
        'std_error': best_model.bse.values,
        'z_value': best_model.tvalues.values,
        'p_value': best_model.pvalues.values,
        'ci_lower': best_model.conf_int()[0].values,
        'ci_upper': best_model.conf_int()[1].values
    })
    coef_table.to_csv(args.coefficients, index=False)
    print(f"Coefficients saved to {args.coefficients}")

    # Pseudo-R²
    pseudo_r2 = 1 - (best_model.deviance / best_model.null_deviance)

    # Save summary
    with open(args.output, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write(f"{model_name} GLM Regression: Introner Effects on Isoform Diversity\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Model: {formula}\n")
        f.write(f"Family: {model_name}\n\n")

        f.write("Data Summary:\n")
        f.write(f"  Sample size: {len(data)}\n")
        f.write(f"  Mean(n_isoforms): {mean_val:.3f}\n")
        f.write(f"  Variance(n_isoforms): {var_val:.3f}\n")
        f.write(f"  VMR: {vmr:.3f}\n\n")

        f.write("Model Fit:\n")
        f.write(f"  AIC: {best_model.aic:.2f}\n")
        f.write(f"  BIC: {best_model.bic:.2f}\n")
        f.write(f"  Deviance: {best_model.deviance:.2f}\n")
        f.write(f"  Pearson chi2: {best_model.pearson_chi2:.2f}\n")
        f.write(f"  Pseudo-R2: {pseudo_r2:.4f}\n\n")

        f.write("Overdispersion Diagnostics:\n")
        f.write(f"  Dispersion (Pearson): {dispersion_pearson:.3f}\n")
        f.write(f"  Dispersion (Deviance): {dispersion_deviance:.3f}\n\n")

        f.write("Full Model Summary:\n")
        f.write("=" * 80 + "\n")
        f.write(str(best_model.summary()))
        f.write("\n\n")

        f.write("Coefficients (Exponentiated):\n")
        f.write("=" * 80 + "\n")
        exp_table = pd.DataFrame({
            'Variable': best_model.params.index,
            'exp(Estimate)': np.exp(best_model.params.values),
            'exp(CI_lower)': np.exp(best_model.conf_int()[0].values),
            'exp(CI_upper)': np.exp(best_model.conf_int()[1].values)
        })
        f.write(exp_table.to_string(index=False))
        f.write("\n\n")

        # Key findings
        f.write("Key Findings:\n")
        f.write("=" * 80 + "\n")
        for var in ['introner_gain', 'introner_loss']:
            if var in best_model.params.index:
                coef = best_model.params[var]
                pval = best_model.pvalues[var]
                exp_coef = np.exp(coef)
                sig = "SIGNIFICANT" if pval < 0.05 else "Not significant"
                f.write(f"{var}: exp(coef)={exp_coef:.4f}, p={pval:.4g} ({sig})\n")

    print(f"Summary saved to {args.output}")

    # Diagnostic plots
    pearson_resid = best_model.resid_pearson
    deviance_resid = best_model.resid_deviance
    fitted_values = best_model.fittedvalues

    with PdfPages(args.diagnostics) as pdf:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        # 1. Deviance Residuals vs Fitted
        ax = axes[0, 0]
        ax.scatter(fitted_values, deviance_resid, alpha=0.3, s=10)
        ax.axhline(y=0, color='r', linestyle='--', linewidth=2)
        try:
            lw = lowess(deviance_resid, fitted_values, frac=0.2)
            ax.plot(lw[:, 0], lw[:, 1], 'b-', linewidth=2)
        except Exception:
            pass
        ax.set_xlabel('Fitted values')
        ax.set_ylabel('Deviance residuals')
        ax.set_title('Deviance Residuals vs Fitted')

        # 2. Q-Q plot
        ax = axes[0, 1]
        stats.probplot(deviance_resid, dist="norm", plot=ax)
        ax.set_title('Normal Q-Q (Deviance Residuals)')

        # 3. Scale-Location
        ax = axes[1, 0]
        resid_abs_sqrt = np.sqrt(np.abs(deviance_resid))
        ax.scatter(fitted_values, resid_abs_sqrt, alpha=0.3, s=10)
        ax.set_xlabel('Fitted values')
        ax.set_ylabel('sqrt(|Deviance residuals|)')
        ax.set_title('Scale-Location')

        # 4. Histogram of residuals
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
