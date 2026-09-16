"""Two-strain R2C2 gene-FE gain/loss fit; retains the R2C2 gene universe."""
import argparse
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from prepare_r2c2_expression import fit, normalize


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--gain-coefficients", type=Path, required=True)
    args = parser.parse_args()
    base = args.data_dir
    source = base / 'r2c2/expression_model_data.tsv'
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    original = pd.read_csv(source, sep='\t')
    data = original[original.strain.isin(['CCMP1545', 'RCC1614'])].copy()
    data, norm, n = normalize(data, True)
    data.to_csv(out/'expression_model_data.tsv', sep='\t', index=False)
    norm.to_csv(out/'normalization_factors.tsv', sep='\t', index=False)
    # A within-gene comparison requires both strains after callability filtering.
    complete = data[data.complete_event_callability].copy()
    complete = complete[complete.groupby('gene_id').strain.transform('nunique').eq(2)].copy()
    assert complete.groupby('gene_id').size().eq(8).all()
    assert set(complete.strain) == {'CCMP1545', 'RCC1614'}
    assert not complete.mating_type_gene.any()
    assert not complete.duplicated(['gene_id','replicate']).any()
    primary = fit.fit_gene_fe_ppml(complete, fit.TURNOVER_VARS, 'offset_median_ratio', 'group1_gain_loss_gene_fe_ppml')
    coefficients = fit.coefficient_frame(primary, 'group1_gain_loss_gene_fe_ppml', 'gene-fixed-effect Poisson pseudo-likelihood')
    assert 'strain_RCC1749' not in primary.params
    assert primary.converged
    assert np.isfinite(coefficients[['coefficient','std_error','p_value']]).all(axis=None)
    assert coefficients.std_error.gt(0).all()
    primary.model_data.to_csv(out/'primary_model_data.tsv',sep='\t',index=False)
    coefficients.to_csv(out/'model_coefficients.tsv',sep='\t',index=False)
    focal = coefficients[coefficients.term.isin(fit.TURNOVER_VARS)].copy()
    focal.to_csv(out/'gain_loss_tests.tsv',sep='\t',index=False)
    support = fit.support_row('group1_gain_loss_gene_fe_ppml',primary.model_data)
    pd.DataFrame([support]).to_csv(out/'model_support.tsv',sep='\t',index=False)
    three = pd.read_csv(args.gain_coefficients,sep='\t')
    three = three[three.model.eq('paired_turnover_gene_fe_ppml') & three.term.isin(fit.TURNOVER_VARS)]
    three.merge(focal,on='term',suffixes=('_three_strains','_group1'),validate='one_to_one').to_csv(
        out/'comparison_with_three_strains.tsv',sep='\t',index=False)
    # Check original observations are unchanged; only normalization is recomputed.
    keys=['gene_id','strain','replicate']
    old=original.set_index(keys).loc[data.set_index(keys).index]
    np.testing.assert_array_equal(old.raw_count,data.raw_count)
    np.testing.assert_allclose(data.offset_median_ratio,np.log(data.median_ratio_size_factor))
    assert abs(np.log(norm.median_ratio_size_factor).mean()) < 1e-10
    lines=['CCMP1545–RCC1614-only R2C2 fixed-effect gain/loss model','',
        'Four biological replicates per strain. Both strains required after complete-event-callability filtering.',
        'Same gain/loss and GC predictors; CCMP1545 reference; gene-clustered sandwich SEs.',
        'Median-ratio normalization recalculated for eight replicates; no length offset.',
        'Retains the R2C2 exclusion universe, including common-ID exclusions derived from RCC1749.',
        'No new high-confidence orthology or current-missing-call filter added; this matches the existing gain/loss feature definitions.',
        f'Normalization genes: {n}; primary genes: {support["genes"]}; observations: {support["observations"]}.',
        f'Gene-strain rows with gain: {support["gain_positive_gene_strain_rows"]}; loss: {support["loss_positive_gene_strain_rows"]}.',
        '']
    for r in focal.itertuples():
        lines.append(f'{r.term}: RR {r.rate_ratio:.6f} ({r.rr_ci_lower:.6f}–{r.rr_ci_upper:.6f}); '
                     f'P={r.p_value:.8g}.')
    lines.append('PASS: convergence, paired observations, unchanged counts, and depth-only offsets.')
    (out/'summary.txt').write_text('\n'.join(lines)+'\n')
    (out/'provenance.json').write_text(json.dumps({'input':str(source),'sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'fit_script_sha256':hashlib.sha256(Path(fit.__file__).read_bytes()).hexdigest()},indent=2))
    print('\n'.join(lines))


if __name__=='__main__': main()
