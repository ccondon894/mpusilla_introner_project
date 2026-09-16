"""Three-strain gene-FE expression association for exactly one varying introner locus."""
import argparse
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from prepare_r2c2_expression import ROOT, fit, normalize

SAMPLES = ['CCMP1545','RCC1614','RCC1749']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    base = args.data_dir
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    paths = {'genotypes': ROOT/'results/genotyping/genotype_matrix.final.tsv',
        'events':ROOT/'results/expression/isoform_analysis/turnover/locus_level_reference_comparisons.tsv',
        'mapping':ROOT/'results/expression/isoform_analysis/turnover/locus_gene_map.tsv',
        'mt':base/'excluded_mating_type_gene_ids.tsv'}
    gt = pd.read_csv(paths['genotypes'],sep='\t',low_memory=False)
    calls = gt.pivot(index='ortholog_id',columns='sample',values='presence')[SAMPLES]
    e = pd.read_csv(paths['events'],sep='\t').rename(columns={'common_gene_id':'gene_id'})
    e = e[e.strain.isin(SAMPLES)].copy()
    assert not e.duplicated(['ortholog_id','strain']).any()
    assert e.groupby('ortholog_id').gene_id.nunique().eq(1).all()
    truth = gt[['ortholog_id','sample','presence']].rename(columns={'sample':'strain'})
    verify = e.merge(truth,on=['ortholog_id','strain'],validate='one_to_one')
    assert len(verify)==len(e)
    assert verify.sample_call.eq(verify.presence).all()
    assert e.reference_call.eq(e.ortholog_id.map(calls.CCMP1545)).all()
    loci = e.groupby('ortholog_id').agg(gene_id=('gene_id','first'),
        link_confident=('event_gene_link_confident','all'),n_strains=('strain','nunique')).join(calls)
    loci['complete_calls'] = loci[SAMPLES].isin([1,2]).all(axis=1) & loci.n_strains.eq(3)
    confidence = gt.groupby('ortholog_id').within_group_orthology_confidence.agg(lambda x:x.eq('high').all())
    loci['high_within_group_orthology'] = loci.index.map(confidence)
    loci['variable_locus'] = loci.complete_calls & loci[SAMPLES].nunique(axis=1).gt(1)
    loci['eligible_locus'] = loci.complete_calls & loci.link_confident & loci.high_within_group_orthology
    bad = set(loci.loc[~loci.eligible_locus,'gene_id'])
    mapping = pd.read_csv(paths['mapping'],sep='\t')
    bad.update(mapping.loc[mapping.ortholog_gene_mapping_status.eq('conflicting_common_genes'),'native_gene_id'].dropna())
    bad.update(pd.read_csv(paths['mt'],sep='\t').gene_id)
    features = loci.groupby('gene_id').agg(n_loci=('variable_locus','size'),
        n_variable_loci=('variable_locus','sum'),all_loci_eligible=('eligible_locus','all'))
    features['eligible_gene'] = features.all_loci_eligible & ~features.index.isin(bad)
    features['single_variable_locus'] = features.eligible_gene & features.n_variable_loci.eq(1)
    loci.to_csv(out/'locus_audit.tsv',sep='\t')
    features.to_csv(out/'gene_eligibility.tsv',sep='\t')
    focal_loci = loci[loci.gene_id.isin(features.index[features.single_variable_locus]) & loci.variable_locus]
    assert not focal_loci.gene_id.duplicated().any()
    focal_loci.to_csv(out/'single_variable_loci.tsv',sep='\t')
    # Retain stable controls only as an explicitly separate sensitivity model.
    single = set(features.index[features.single_variable_locus])
    stable = set(features.index[features.eligible_gene & features.n_variable_loci.eq(0)])
    coefficients, supports, used_keys = [], [], {}
    for platform in ['r2c2','short_read_matched']:
        path = base/platform/'expression_model_data.tsv'
        paths[platform] = path
        data = pd.read_csv(path,sep='\t')
        data,norm,n = normalize(data,platform=='r2c2')
        norm.to_csv(out/f'{platform}.normalization.tsv',sep='\t',index=False)
        data = data[data.groupby('gene_id').strain.transform('nunique').eq(3)].copy()
        for subset,genes in [('single_locus_only',single),('single_locus_plus_stable',single|stable)]:
            frame = data[data.gene_id.isin(genes)].copy()
            # Reconstruct burden from verified locus calls, rather than trusting aggregate counts.
            for strain in SAMPLES:
                burden = loci.groupby('gene_id')[strain].agg(lambda x:int(x.eq(1).sum()))
                rows=frame.strain.eq(strain)
                np.testing.assert_array_equal(frame.loc[rows,'current_introner_count'],frame.loc[rows,'gene_id'].map(burden))
            frame['focal_present'] = 0
            for strain in SAMPLES:
                rows=frame.strain.eq(strain)
                assignment=focal_loci.set_index('gene_id')[strain].eq(1).astype(int)
                frame.loc[rows,'focal_present']=frame.loc[rows,'gene_id'].map(assignment).fillna(0).astype(int)
            assert frame.groupby('gene_id').size().eq(12).all()
            assert not frame.mating_type_gene.any()
            assert frame.current_missing_locus_count.eq(0).all()
            model=f'{platform}_{subset}'
            # Explicitly ensure within-gene strain, presence and GC design is identifiable.
            design=pd.DataFrame({'strain_RCC1614':frame.strain.eq('RCC1614').astype(float),
                'strain_RCC1749':frame.strain.eq('RCC1749').astype(float),
                'focal_present':frame.focal_present,'GC_content':frame.GC_content})
            centered=design-design.groupby(frame.gene_id).transform('mean')
            assert np.linalg.matrix_rank(centered.to_numpy())==4, 'Nonidentifiable single-locus design'
            result=fit.fit_gene_fe_ppml(frame,['focal_present'],'offset_median_ratio',model)
            used=result.model_data
            keys=used[['gene_id','strain','replicate']].sort_values(['gene_id','strain','replicate']).reset_index(drop=True)
            if platform=='r2c2':used_keys[subset]=keys
            else:pd.testing.assert_frame_equal(keys,used_keys[subset])
            assert result.converged
            coef=fit.coefficient_frame(result,model,'gene-fixed-effect Poisson pseudo-likelihood')
            assert np.isfinite(coef[['coefficient','std_error','p_value']]).all(axis=None)
            assert coef.std_error.gt(0).all()
            coef['platform']=platform;coef['subset']=subset
            coefficients.append(coef)
            used.to_csv(out/f'{model}.data.tsv',sep='\t',index=False)
            supports.append(dict(model=model,genes=used.gene_id.nunique(),observations=len(used),
                single_locus_genes=used.loc[used.gene_id.isin(single),'gene_id'].nunique(),converged=result.converged))
    coefficients=pd.concat(coefficients,ignore_index=True)
    coefficients.to_csv(out/'model_coefficients.tsv',sep='\t',index=False)
    effects=coefficients[coefficients.term.eq('focal_present')].copy()
    effects.to_csv(out/'presence_effects.tsv',sep='\t',index=False)
    pd.DataFrame(supports).to_csv(out/'model_support.tsv',sep='\t',index=False)
    used_single=set(pd.read_csv(out/'r2c2_single_locus_only.data.tsv',sep='\t').gene_id)
    patterns=focal_loci[focal_loci.gene_id.isin(used_single)][SAMPLES].replace({1:'P',2:'A'}).agg('/'.join,axis=1).value_counts()
    patterns.rename_axis('CCMP1545_RCC1614_RCC1749').to_csv(out/'presence_patterns.tsv',sep='\t',header=['genes'])
    lines=['Single variable introner locus per gene across CCMP1545, RCC1614 and RCC1749','',
        'Primary: exactly one variable locus across all three strains; every other mapped locus callable and invariant.',
        'All loci: high within-group orthology confidence and confident event-to-gene linkage; conflicting candidate genes excluded.',
        'Cross-group origin/confidence labels describe ancestral/independent origin, not a newly validated cross-group orthology filter.',
        'Calls checked against final genotype matrix; absent=2 regardless of interval length; missing=3 never treated as absent.',
        'All three strains and four expression replicates per strain required; MT exclusions retained.',
        'Gene-FE PPML: strain + focal_present + GC + log(exposure), gene-clustered uncertainty.',
        'Focal presence RR compares carrying vs lacking the one variable introner; shared introners absorbed by gene FE.',
        'Short-read depth/exon-length exposure; R2C2 depth-only exposure. Same observation keys across platforms.',
        'Exploratory, not independent replication; mapped-locus completeness cannot exclude undiscovered/unmapped introners.',
        '',pd.DataFrame(supports).to_string(index=False),'',patterns.to_string(),'']
    for r in effects.itertuples():lines.append(f'{r.model}: RR {r.rate_ratio:.5f} ({r.rr_ci_lower:.5f}–{r.rr_ci_upper:.5f}); P={r.p_value:.6g}')
    lines.append('PASS: genotype concordance, reconstructed burden, full-rank design, matched observations, convergence and finite estimates.')
    (out/'summary.txt').write_text('\n'.join(lines)+'\n')
    paths['script']=Path(__file__);paths['fitter']=Path(fit.__file__)
    (out/'provenance.json').write_text(json.dumps({k:{'path':str(v),'sha256':hashlib.sha256(v.read_bytes()).hexdigest()} for k,v in paths.items()},indent=2))
    print('\n'.join(lines))


if __name__=='__main__':main()
