#!/usr/bin/env python3
"""Fit the R2C2 within/between-gene NB-GEE association reported in Figure 5B."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from prepare_r2c2_expression import fit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    data = pd.read_csv(args.input, sep='\t')
    data = fit.add_within_between(data, fit.TURNOVER_VARS + fit.STATE_VARS + ['GC_content'])
    paired = data[data.groupby('gene_id').strain.transform('nunique').ge(2)].copy()
    complete = paired[paired.complete_event_callability].copy()
    alpha = fit.estimate_within_replicate_alpha(complete, 'offset_median_ratio')
    formula = ("raw_count ~ C(strain, Treatment(reference='CCMP1545'))"
               " + current_introner_count_within + current_introner_count_between"
               " + GC_content_within + GC_content_between")
    name = 'current_state_nb_gee'
    result = fit.fit_gee(paired, formula, 'offset_median_ratio', alpha, name)
    coefficients = fit.coefficient_frame(result, name, 'negative binomial GEE')
    if not result.converged or not np.isfinite(coefficients[['coefficient','std_error','p_value','rate_ratio']]).all(axis=None):
        raise RuntimeError('R2C2 NB-GEE model failed validation')
    args.outdir.mkdir(parents=True, exist_ok=True)
    coefficients.to_csv(args.outdir/'model_coefficients.tsv',sep='\t',index=False)
    pd.DataFrame([fit.support_row(name, paired, alpha, float(result.cov_struct.dep_params))]).to_csv(args.outdir/'model_support.tsv',sep='\t',index=False)
    focal = coefficients[coefficients.term.eq('current_introner_count_between')].iloc[0]
    (args.outdir/'summary.txt').write_text(
        f'Between-gene R2C2 association\nRR = {focal.rate_ratio:.12g}; '
        f'95% CI {focal.rr_ci_lower:.12g}–{focal.rr_ci_upper:.12g}; P = {focal.p_value:.12g}\n')


if __name__ == '__main__':
    main()
