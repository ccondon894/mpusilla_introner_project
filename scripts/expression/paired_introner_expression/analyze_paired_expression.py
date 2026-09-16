#!/usr/bin/env python3
"""CCMP1545/RCC1614 gene-paired introner presence analysis (no ancestral polarity)."""
from pathlib import Path
import json
import hashlib
import os
os.environ.setdefault('MPLCONFIGDIR', '/scratch1/chris/tmp/matplotlib')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'results/expression/functional/paired_introner_expression'
REF, ALT = 'CCMP1545', 'RCC1614'


def fit_pair(table, adjusted=True):
    """Profile gene intercepts from PPML; independent gene sandwich covariance.

    Conditional on each gene's total, the profiled mean allocates that total
    between strains using a logistic probability. This is a Poisson mean-model
    estimator, not an assumption that reads are independent binomial trials.
    """
    terms = ['strain_RCC1614', 'delta_present_count']
    x = np.column_stack([np.ones(len(table)), table.delta_present_count])
    if adjusted:
        terms.append('delta_GC')
        x = np.column_stack([x, table.delta_GC])
    if np.linalg.matrix_rank(x) != x.shape[1]:
        raise ValueError('Nonidentifiable paired design')
    n = table.total_count.to_numpy(float)
    y = table.count_alt.to_numpy(float)
    off = table.log_exposure_ratio.to_numpy(float)
    scale = n.sum()
    def objective(b):
        eta = off + x @ b
        return (np.sum(n*np.logaddexp(0, eta)-y*eta)/scale,
                x.T @ (n*expit(eta)-y)/scale)
    opt = minimize(objective, np.zeros(x.shape[1]), jac=True, method='BFGS',
                   options={'gtol': 1e-11, 'maxiter': 1000})
    assert np.max(np.abs(opt.jac)) < 1e-7, opt.message
    p = expit(off+x@opt.x)
    bread = x.T @ ((n*p*(1-p))[:,None]*x)
    score = x*(y-n*p)[:,None]
    inv = np.linalg.inv(bread)
    covariance = inv @ (score.T@score) @ inv * len(x)/(len(x)-x.shape[1])
    se = np.sqrt(np.diag(covariance))
    result = pd.DataFrame({'term':terms, 'coefficient':opt.x, 'std_error':se})
    result['p_value'] = 2*norm.sf(np.abs(opt.x/se))
    result['ci_lower'] = opt.x-1.96*se
    result['ci_upper'] = opt.x+1.96*se
    for target, source in [('rate_ratio','coefficient'), ('rr_lower','ci_lower'), ('rr_upper','ci_upper')]:
        result[target] = np.exp(result[source])
    return result


def validate_estimator():
    """Check the profiled estimator against explicit gene-intercept GLM."""
    rng = np.random.default_rng(923)
    k = 60
    delta = rng.integers(-2, 3, k)
    gc = rng.normal(0, .02, k)
    exposure = rng.uniform(.7, 1.5, (k,2))
    baseline = rng.uniform(20, 200, k)
    count0 = rng.poisson(baseline*exposure[:,0])
    count1 = rng.poisson(baseline*exposure[:,1]*np.exp(.2-.3*delta+gc))
    t = pd.DataFrame({'delta_present_count':delta, 'delta_GC':gc,
        'total_count':count0+count1, 'count_alt':count1,
        'log_exposure_ratio':np.log(exposure[:,1]/exposure[:,0])})
    got = fit_pair(t).coefficient.to_numpy()
    design = np.column_stack([np.repeat(np.eye(k),2,axis=0),
        np.tile([0,1],k), np.ravel(np.column_stack([np.zeros(k),delta])),
        np.ravel(np.column_stack([np.zeros(k),gc]))])
    expected = sm.GLM(np.column_stack([count0,count1]).ravel(),design,
        offset=np.log(exposure).ravel(),family=sm.families.Poisson()).fit().params[-3:]
    assert np.allclose(got,expected,atol=2e-5), (got,expected)
    reverse = t.copy()
    reverse['count_alt'] = t.total_count-t.count_alt
    reverse['delta_present_count'] *= -1
    reverse['delta_GC'] *= -1
    reverse['log_exposure_ratio'] *= -1
    rev = fit_pair(reverse).coefficient.to_numpy()
    assert np.allclose(got[1:],rev[1:],atol=2e-5)
    assert np.isclose(got[0],-rev[0],atol=2e-5)


def main():
    global OUT, ROOT
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=ROOT)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    ROOT, OUT = args.project_root.resolve(), args.outdir
    OUT.mkdir(parents=True, exist_ok=True)
    validate_estimator()
    inputs = {
        'counts': ROOT/'results/expression/functional/short_read_data/expression_model_data.tsv',
        'genotypes': ROOT/'results/genotyping/genotype_matrix.final.tsv',
        'events': ROOT/'results/expression/isoform_analysis/turnover/locus_level_reference_comparisons.tsv',
        'mapping': ROOT/'results/expression/isoform_analysis/turnover/locus_gene_map.tsv',
        'gtf': ROOT/'results/annotations/CCMP1545.gtf',
    }
    counts = pd.read_csv(inputs['counts'],sep='\t')
    counts = counts[counts.strain.isin([REF,ALT])].copy()
    assert not counts.duplicated(['gene_id','replicate']).any()
    assert counts.groupby('strain').replicate.nunique().eq(4).all()
    events = pd.read_csv(inputs['events'],sep='\t')
    events = events[events.strain.eq(ALT)].rename(columns={'common_gene_id':'gene_id'})
    assert not events.ortholog_id.duplicated().any()
    gt = pd.read_csv(inputs['genotypes'],sep='\t',low_memory=False)
    calls = gt.pivot(index='ortholog_id',columns='sample',values='presence')
    events = events.merge(calls[[REF,ALT]],left_on='ortholog_id',right_index=True,validate='one_to_one')
    assert events.reference_call.eq(events[REF]).all()
    assert events.sample_call.eq(events[ALT]).all()
    confidence = gt.groupby('ortholog_id').within_group_orthology_confidence.agg(lambda x: (x=='high').all())
    events['high_orthology'] = events.ortholog_id.map(confidence)
    events['callable_pair'] = events[REF].isin([1,2]) & events[ALT].isin([1,2])
    events['discordant'] = events.callable_pair & events[REF].ne(events[ALT])
    events['eligible_locus'] = (events.callable_pair & events.high_orthology &
        events.event_gene_link_confident & ~events.mating_type_gene)
    events.to_csv(OUT/'locus_audit.tsv',sep='\t',index=False)
    bad_genes = set(events.loc[~events.eligible_locus,'gene_id'])
    mapping = pd.read_csv(inputs['mapping'],sep='\t')
    # Ambiguous loci can otherwise make their candidate genes look introner-free.
    bad_genes.update(mapping.loc[mapping.ortholog_gene_mapping_status.eq('conflicting_common_genes'),
                                 'native_gene_id'].dropna())
    # Independently exclude all reference genes overlapping the mating-type region.
    import sys
    sys.path.insert(0,str(ROOT/'scripts/expression/negative_binomial_expression_pilot'))
    from prepare_expression_pilot_data import reference_gene_table
    bounds = reference_gene_table(inputs['gtf'])
    mt = set(bounds.loc[bounds.contig.eq('CCMP1545#0#scaffold_2') &
        bounds.gene_start.le(1730591) & bounds.gene_end.ge(49808),'gene_id'])
    bad_genes.update(mt)
    features = events.groupby('gene_id').agg(
        n_loci=('ortholog_id','size'), n_discordant=('discordant','sum'),
        present_ref=(REF,lambda x:int(x.eq(1).sum())),
        present_alt=(ALT,lambda x:int(x.eq(1).sum())))
    counts['exposure'] = counts.median_ratio_size_factor * counts.effective_exon_length/1000
    assert np.isfinite(counts.exposure).all() and counts.exposure.gt(0).all()
    aggregate = counts.groupby(['gene_id','strain']).agg(count=('raw_count','sum'),
        exposure=('exposure','sum'), gc=('GC_content','first'),
        length=('effective_exon_length','first'), n_replicates=('replicate','nunique'))
    wide = aggregate.unstack('strain').dropna()
    pair = pd.DataFrame(index=wide.index)
    for name in ['count','exposure','gc','length','n_replicates']:
        pair[name+'_ref'] = wide[name][REF]
        pair[name+'_alt'] = wide[name][ALT]
    pair = pair.join(features).fillna({'n_loci':0,'n_discordant':0,'present_ref':0,'present_alt':0})
    pair['delta_present_count'] = pair.present_alt-pair.present_ref
    pair['delta_GC'] = pair.gc_alt-pair.gc_ref
    pair['log_exposure_ratio'] = np.log(pair.exposure_alt/pair.exposure_ref)
    pair['total_count'] = pair.count_alt+pair.count_ref
    pair['excluded_callability_linkage_or_MT'] = pair.index.isin(bad_genes)
    pair['eligible'] = (~pair.excluded_callability_linkage_or_MT & pair.total_count.gt(0)
        & pair.n_replicates_ref.eq(4) & pair.n_replicates_alt.eq(4))
    pair['log_expression_ratio_pseudocount'] = np.log((pair.count_alt+.5)/pair.exposure_alt)-np.log((pair.count_ref+.5)/pair.exposure_ref)
    pair.to_csv(OUT/'gene_pair_audit.tsv',sep='\t')
    primary = pair[pair.eligible].copy()
    assert not primary.index.isin(mt).any()
    assert primary.n_discordant.ge(primary.delta_present_count.abs()).all()
    discordant = calls[calls[REF].isin([1,2]) & calls[ALT].isin([1,2]) & calls[REF].ne(calls[ALT])][[REF,ALT]].copy()
    discordant = discordant.join(events.set_index('ortholog_id')[['gene_id','eligible_locus']])
    discordant['included'] = discordant.gene_id.isin(primary.index)
    discordant['exclusion_reason'] = np.select(
        [discordant.included,discordant.gene_id.isna(),~discordant.gene_id.isin(pair.index),
         discordant.gene_id.isin(bad_genes)],
        ['included','no_unique_common_gene_in_event_table','no_complete_expression_pair',
         'gene_excluded_for_callability_linkage_orthology_or_MT'],
        default='zero_total_counts_or_incomplete_replicates')
    assert int(discordant.included.sum()) == int(primary.n_discordant.sum())
    discordant.to_csv(OUT/'all_discordant_locus_inclusion.tsv',sep='\t')
    subsets = {
        'primary_all_callable_genes': (primary, True),
        'single_discordant_plus_stable_controls': (primary[primary.n_discordant.le(1)], True),
        'discordant_genes_only': (primary[primary.n_discordant.gt(0)], True),
        'single_discordant_genes_only': (primary[primary.n_discordant.eq(1)], True),
        'without_GC_adjustment': (primary, False),
        'exonic_length_within_10_percent': (primary[(primary.length_alt/primary.length_ref).between(1/1.1,1.1)], True),
    }
    coefficients, support = [], []
    for name,(data,adjusted) in subsets.items():
        result = fit_pair(data,adjusted)
        result.insert(0,'model',name)
        coefficients.append(result)
        support.append({'model':name,'genes':len(data),
            'discordant_genes':int(data.n_discordant.gt(0).sum()),
            'discordant_loci':int(data.n_discordant.sum()),
            'more_present_in_ref':int(data.delta_present_count.lt(0).sum()),
            'more_present_in_alt':int(data.delta_present_count.gt(0).sum()),
            'equal_burden':int(data.delta_present_count.eq(0).sum())})
    coef = pd.concat(coefficients,ignore_index=True)
    coef.to_csv(OUT/'model_coefficients.tsv',sep='\t',index=False)
    pd.DataFrame(support).to_csv(OUT/'model_support.tsv',sep='\t',index=False)
    # Equal gene weighting checks dependence on the count-weighted PPML estimand.
    x = sm.add_constant(primary[['delta_present_count','delta_GC']])
    logfit = sm.OLS(primary.log_expression_ratio_pseudocount,x).fit(cov_type='HC3')
    pd.DataFrame({'term':logfit.params.index,'coefficient':logfit.params.values,
        'std_error':logfit.bse.values,'p_value':logfit.pvalues.values,
        'ci_lower':logfit.conf_int()[0].values,'ci_upper':logfit.conf_int()[1].values}).to_csv(
            OUT/'equal_gene_log_ratio_sensitivity.tsv',sep='\t',index=False)
    chosen = coef[(coef.model=='primary_all_callable_genes')].set_index('term')
    primary['adjusted_log_ratio_alt_ref'] = primary.log_expression_ratio_pseudocount - chosen.loc['strain_RCC1614','coefficient'] - primary.delta_GC*chosen.loc['delta_GC','coefficient']
    primary['present_vs_absent_log_ratio_single_locus'] = np.where(primary.n_discordant.eq(1),
        primary.adjusted_log_ratio_alt_ref*np.sign(primary.delta_present_count),np.nan)
    primary.to_csv(OUT/'paired_expression_results.tsv',sep='\t')
    # Forest plot: all prespecified comparisons, including null/opposite estimates.
    effects = coef[coef.term.eq('delta_present_count')].copy()
    labels = ['Primary: all callable genes','Single difference + stable controls',
        'Discordant genes only','Single difference genes only','Without GC adjustment',
        'Exonic length differs by ≤10%']
    fig,ax = plt.subplots(figsize=(8,4.5))
    y = np.arange(len(effects))
    ax.errorbar(effects.rate_ratio,y,xerr=np.vstack([effects.rate_ratio-effects.rr_lower,
        effects.rr_upper-effects.rate_ratio]),fmt='o',color='#0072b2',capsize=3)
    ax.axvline(1,color='#999999',ls='--')
    ax.set_yticks(y,labels);ax.invert_yaxis();ax.set_xscale('log')
    ax.set_xlabel('Expression ratio per additional present introner (95% CI)')
    ax.set_title('CCMP1545–RCC1614: gene-paired expression')
    fig.tight_layout();fig.savefig(OUT/'paired_expression_effects.png',dpi=220);plt.close(fig)
    manifest = {k:{'path':str(v),'sha256':hashlib.sha256(v.read_bytes()).hexdigest()} for k,v in inputs.items()}
    (OUT/'input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (OUT/'validation_summary.txt').write_text(
        'PASS: profile estimator matches explicit gene fixed effects on synthetic counts\n'
        'PASS: reversing strain labels preserves introner coefficient\n'
        'PASS: event calls match current final genotype matrix\n'
        'PASS: four replicates per strain and unique gene/replicate keys\n'
        'PASS: positive finite exposures and full-rank converged models\n'
        'PASS: reference mating-type-overlapping genes independently excluded\n'
        'PASS: missing/uncertain mapped loci exclude whole genes\n')
    lines = ['CCMP1545–RCC1614 paired expression analysis', '',
        f'Eligible genes: {len(primary)}; discordant genes: {int(primary.n_discordant.gt(0).sum())}; discordant loci: {int(primary.n_discordant.sum())}', '']
    for r in effects.itertuples():
        lines.append(f'{r.model}: RR {r.rate_ratio:.4f} (95% CI {r.rr_lower:.4f}–{r.rr_upper:.4f}), P={r.p_value:.5g}')
    lines += ['',f'Equal-gene log-ratio sensitivity: beta={logfit.params["delta_present_count"]:.4f}, P={logfit.pvalues["delta_present_count"]:.5g}',
        '', 'All comparisons are exploratory. No ancestral direction is inferred. No causal interpretation.',
        'Confidence intervals use independent genes, not independent strains; population-wide generalization is limited.']
    (OUT/'results_summary.txt').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
