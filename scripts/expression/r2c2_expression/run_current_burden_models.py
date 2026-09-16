"""Within-gene burden association using R2C2 or matched short-read counts."""
import argparse
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from prepare_r2c2_expression import fit, normalize


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--count-source', choices=['r2c2','short_read_matched'], default='r2c2')
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    long_reads = args.count_source == 'r2c2'
    base = args.data_dir
    source = base/args.count_source/'expression_model_data.tsv'
    original = pd.read_csv(source, sep='\t')
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    coefficients, support, checks = [], [], []
    for name, strains in [('three_strains',['CCMP1545','RCC1614','RCC1749']),
                          ('group1_only',['CCMP1545','RCC1614'])]:
        data = original[original.strain.isin(strains)].copy()
        data, normalization, n = normalize(data, long_reads)
        normalization.to_csv(out/f'{name}.normalization.tsv',sep='\t',index=False)
        assert not data.mating_type_gene.any()
        assert not data.duplicated(['gene_id','replicate']).any()
        length = 1.0 if long_reads else data.effective_exon_length / 1000
        np.testing.assert_allclose(data.offset_median_ratio,np.log(data.median_ratio_size_factor * length))
        assert abs(np.log(normalization.median_ratio_size_factor).mean()) < 1e-10
        for eligibility in ['complete_current_calls','previous_eligibility']:
            frame = data.copy()
            if eligibility == 'complete_current_calls':
                frame = frame[frame.current_missing_locus_count.eq(0)].copy()
            frame = frame[frame.groupby('gene_id').strain.transform('nunique').ge(2)].copy()
            model = f'{name}_{eligibility}'
            result = fit.fit_gene_fe_ppml(frame,['current_introner_count'],'offset_median_ratio',model)
            c = fit.coefficient_frame(result,model,'gene-fixed-effect Poisson pseudo-likelihood')
            c['strain_set'] = name
            c['eligibility'] = eligibility
            coefficients.append(c)
            used = result.model_data
            used.to_csv(out/f'{model}.data.tsv',sep='\t',index=False)
            varying = used.groupby('gene_id').current_introner_count.nunique().gt(1).sum()
            support.append(dict(model=model,genes=used.gene_id.nunique(),observations=len(used),
                varying_burden_genes=int(varying),normalization_genes=n,
                missing_current_gene_strain_rows=int(used.drop_duplicates(['gene_id','strain']).current_missing_locus_count.gt(0).sum()),
                converged=result.converged))
            assert result.converged
            assert np.isfinite(c[['coefficient','std_error','p_value']]).all(axis=None)
            assert c.std_error.gt(0).all()
            assert used.groupby(['gene_id','strain']).size().eq(4).all()
            if name == 'group1_only':
                assert 'strain_RCC1749' not in result.params
                assert used.groupby('gene_id').size().eq(8).all()
            if eligibility == 'complete_current_calls': assert used.current_missing_locus_count.eq(0).all()
    coef = pd.concat(coefficients,ignore_index=True)
    coef.to_csv(out/'model_coefficients.tsv',sep='\t',index=False)
    pd.DataFrame(support).to_csv(out/'model_support.tsv',sep='\t',index=False)
    focal = coef[coef.term.eq('current_introner_count')].copy()
    focal.to_csv(out/'burden_effects.tsv',sep='\t',index=False)
    checks += ['PASS: all fits converged; finite coefficients and positive standard errors',
               'PASS: four biological replicates per retained gene-strain; paired gene support',
               'PASS: no current missing calls in primary fits; MT flags absent',
               'PASS: count-source-appropriate offsets and geometric-mean-one size factors']
    (out/'validation_summary.txt').write_text('\n'.join(checks)+'\n')
    lines = [f'Within-gene current introner count: {args.count_source}','',
        'Mean model: gene FE + strain + current_introner_count + GC + log(exposure).',
        'Exposure: median-ratio size factor' + ('.' if long_reads else ' times effective exon length in kb.'),
        'Gene-clustered sandwich uncertainty; two-sided tests; no isoform covariate.',
        'Primary: current_missing_locus_count=0, then require >=2 strains per gene.',
        'Sensitivity: previous current-burden eligibility, including incomplete current calls.',
        'Same R2C2 gene universe; normalization recalculated for each strain set.',
        'Genes lacking mapped introner loci retain inherited zero burden; no new locus/orthology audit.',
        ]
    for r in focal.itertuples():
        lines.append(f'{r.model}: RR {r.rate_ratio:.6f} ({r.rr_ci_lower:.6f}–{r.rr_ci_upper:.6f}); '
                     f'P={r.p_value:.8g}.')
    lines += ['',pd.DataFrame(support).to_string(index=False)]
    (out/'summary.txt').write_text('\n'.join(lines)+'\n')
    (out/'provenance.json').write_text(json.dumps({
        'input':str(source),'sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'fitter_sha256':hashlib.sha256(Path(fit.__file__).read_bytes()).hexdigest()},indent=2))
    print('\n'.join(lines))


if __name__ == '__main__': main()
