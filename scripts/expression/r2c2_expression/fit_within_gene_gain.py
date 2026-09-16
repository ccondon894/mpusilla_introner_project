#!/usr/bin/env python3
"""Fit the within-gene R2C2 gain effect with gene-clustered uncertainty."""
from pathlib import Path
import argparse
import sys
import numpy as np
import pandas as pd
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from figure_color_guide import load_color_guide
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'negative_binomial_expression_pilot'))
from fit_expression_pilot_models import TURNOVER_VARS, fit_gene_fe_ppml, coefficient_frame, support_row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    colors = load_color_guide()
    data = pd.read_csv(args.input, sep='\t')
    paired = data[data.groupby('gene_id').strain.transform('nunique').ge(2)].copy()
    complete = paired[paired.complete_event_callability].copy()
    name = 'paired_turnover_gene_fe_ppml'
    result = fit_gene_fe_ppml(complete, TURNOVER_VARS, 'offset_median_ratio', name)
    coefficients = coefficient_frame(result, name, 'gene-fixed-effect Poisson pseudo-likelihood')
    if not result.converged or not np.isfinite(coefficients[['coefficient','std_error','p_value']]).all(axis=None):
        raise RuntimeError('Within-gene gain fit failed validation')
    gain = coefficients[coefficients.term.eq('reference_relative_gain_count')]
    args.outdir.mkdir(parents=True, exist_ok=True)
    coefficients.to_csv(args.outdir/'model_coefficients.tsv', sep='\t', index=False)
    gain.to_csv(args.outdir/'gain_effect.tsv', sep='\t', index=False)
    result.model_data.to_csv(args.outdir/'model_data.tsv', sep='\t', index=False)
    pd.DataFrame([support_row(name, result.model_data)]).to_csv(args.outdir/'model_support.tsv', sep='\t', index=False)
    r = gain.iloc[0]
    (args.outdir/'summary.txt').write_text(
        f'Within-gene R2C2 gain effect\nRR = {r.rate_ratio:.12g}; 95% CI {r.rr_ci_lower:.12g}–{r.rr_ci_upper:.12g}; P = {r.p_value:.12g}\n'
        f'Gene fixed effects; gain, loss, strain and GC predictors; median-ratio size-factor offset.\n'
        f'Gene-clustered sandwich standard errors; converged in {result.iterations} iterations.\n')
    fig, ax = plt.subplots(figsize=(5, 2.4))
    ax.errorbar(r.rate_ratio, 0, xerr=[[r.rate_ratio-r.rr_ci_lower],[r.rr_ci_upper-r.rate_ratio]], fmt='o', capsize=4, color=colors['Population 1'])
    ax.axvline(1, color=colors['Ancestor'], linestyle='--'); ax.set_xscale('log')
    ax.set_yticks([0], ['Gain, within gene']); ax.set_xlabel('Expression rate ratio (95% CI)')
    fig.tight_layout(); fig.savefig(args.outdir/'within_gene_gain.png', dpi=300); plt.close(fig)


if __name__ == '__main__':
    main()
