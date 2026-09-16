#!/usr/bin/env python3
"""Re-extract all ancestral candidate bodies regionally, reapply QC, and compare legacy estimates."""
from pathlib import Path
import argparse
import concurrent.futures as futures
import importlib.util
import json
import subprocess
import tempfile

import pandas as pd
import pysam
from Bio import SeqIO
from Bio.Seq import Seq

import run_analysis as a


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parent/'results/introner_regional_sensitivity')
    out=parser.parse_args().output
    out.mkdir(parents=True, exist_ok=True)
    (out/'alignments').mkdir(exist_ok=True)
    old = pd.read_csv(a.ROOT/'results/evolution/diversity_metrics/shared_introner_body_dxy_200bp.tsv',sep='\t')
    old = old[old.ancestry_class=='ancestral']
    matrix = pd.read_csv(a.ROOT/'results/genotyping/genotype_matrix.final.tsv',sep='\t')
    matrix = matrix[matrix.ortholog_id.isin(old.ortholog_id)]
    qc = pd.read_csv(a.ROOT/'results/evolution/introner_body_consensus/introner_body_consensus.qc.tsv',sep='\t')
    qc = qc[qc.ortholog_id.isin(old.ortholog_id)].copy()
    spec = importlib.util.spec_from_file_location('qc',a.ROOT/'scripts/evolution/qc_introner_body_consensus.py')
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    seqs = {}
    with tempfile.TemporaryDirectory(prefix='introner_regional_check_',dir=a.TMP) as td, (out/'commands.log').open('w') as log:
        td = Path(td)
        for s in a.SAMPLES:
            rows = matrix[matrix['sample']==s].sort_values(['contig','start'])
            fa = td/f'{s}.fa'
            with pysam.FastaFile(str(a.assembly(s))) as ref, fa.open('w') as h:
                for r in rows.itertuples():
                    start,end=int(r.start)+100,int(r.end)-100
                    h.write(f'>{r.contig}:{start+1}-{end}\n{ref.fetch(r.contig,start,end)}\n')
            result = fa
            if s!='CCMP1545':
                result = td/f'{s}.cons.fa'
                a.run(['bcftools','consensus','-a','N','-f',fa,'-o',result,
                       a.ROOT/f'results/evolution/introner_body_consensus/{s}.introner_body.filtered.vcf.gz'],log)
            with result.open() as handle:
                regional={r.id:str(r.seq).upper() for r in SeqIO.parse(handle,'fasta')}
            with (out/f'{s}.consensus.fa').open('w') as h:
                for r in rows.itertuples():
                    start,end=int(r.start)+100,int(r.end)-100
                    seq=regional[f'{r.contig}:{start+1}-{end}']
                    ix=qc.index[(qc['sample']==s)&(qc.ortholog_id==r.ortholog_id)]
                    assert len(ix)==1
                    ix=ix[0]
                    qc.loc[ix,'consensus_len']=len(seq)
                    qc.loc[ix,'length_ok']=len(seq)==end-start
                    qc.loc[ix,'n_fraction']=seq.count('N')/len(seq) if seq else 1.
                    qc.loc[ix,'consensus_missing']=False
                    if r.orientation=='reverse':
                        seq=str(Seq(seq).reverse_complement())
                    seqs.setdefault(r.ortholog_id,{})[s]=seq
                    h.write(f'>{r.ortholog_id}\n{seq}\n')
    qc['qc_fail_reason']=[mod.fail_reasons(r) for r in qc.to_dict('records')]
    qc['qc_pass']=qc.qc_fail_reason.eq('PASS')
    qc.to_csv(out/'regional_consensus_qc.tsv',sep='\t',index=False)
    counts=qc[qc.qc_pass].groupby('ortholog_id')['sample'].nunique()
    keep=set(counts[counts==13].index)
    results=[]
    with futures.ThreadPoolExecutor(max_workers=2) as pool:
        jobs=[pool.submit(a.align_canonical,lid,seqs[lid],out) for lid in sorted(keep)]
        for job in futures.as_completed(jobs):
            lid,metrics,pairs=job.result()
            results.append(dict(ortholog_id=lid,**metrics))
    df=old[['ortholog_id','dxy_introner']].merge(pd.DataFrame(results),on='ortholog_id',how='left')
    df['regional_minus_legacy']=df.dxy-df.dxy_introner
    df.to_csv(out/'legacy_vs_regional_dxy.tsv',sep='\t',index=False)
    summary=dict(n_candidates=len(old),n_legacy=int(old.dxy_introner.notna().sum()),n_regional=len(keep),
                 n_rescued=int((df.dxy_introner.isna() & df.dxy.notna()).sum()),
                 n_legacy_dropped=int((df.dxy_introner.notna() & df.dxy.isna()).sum()),
                 legacy_mean=float(old.dxy_introner.mean()),
                 regional_mean=float(df.dxy.mean()),mean_change_matched_loci=float(df.regional_minus_legacy.mean()),
                 max_absolute_change=float(df.regional_minus_legacy.abs().max()))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    provenance=dict(script=a.fingerprint(__file__,True),
                    filtered_vcfs=[a.fingerprint(a.ROOT/f'results/evolution/introner_body_consensus/{s}.introner_body.filtered.vcf.gz',True) for s in a.SAMPLES if s!='CCMP1545'])
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    main()
