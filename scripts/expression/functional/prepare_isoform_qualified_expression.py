#!/usr/bin/env python3
"""Prepare CPM summaries for the isoform-qualified expression-class comparisons."""
import argparse
import re
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sqanti', type=Path, required=True)
    p.add_argument('--counts', type=Path, required=True)
    p.add_argument('--gtf', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    sqanti = pd.read_csv(args.sqanti, sep='\t')
    categories = ['full-splice_match', 'novel_in_catalog', 'novel_not_in_catalog']
    richness = sqanti[sqanti.structural_category.isin(categories)].groupby(['gene_id','strain']).isoform.nunique().rename('n_isoforms').reset_index()
    counts = pd.read_csv(args.counts, index_col=0)
    reps = [c for c in counts if c.startswith(('834_', '1614_', '1749_'))]
    cpm = counts[reps].div(counts[reps].sum(), axis=1)*1e6
    cds, mt = defaultdict(int), set()
    for line in args.gtf.open():
        if line.startswith('#'): continue
        fields = line.rstrip().split('\t')
        if len(fields)!=9: continue
        match = re.search(r'gene_id[=\s]+"?([^;"]+)"?', fields[8])
        if not match: continue
        gene = match.group(1)
        if fields[2]=='CDS': cds[gene] += int(fields[4])-int(fields[3])+1
        if fields[0]=='CCMP1545#0#scaffold_2' and int(fields[3])<=1730591 and int(fields[4])>=49808: mt.add(gene)
    frames=[]
    for prefix,strain in [('834_', 'CCMP1545'),('1614_', 'RCC1614'),('1749_', 'RCC1749')]:
        values = cpm[[c for c in reps if c.startswith(prefix)]].mean(axis=1)
        frames.append(pd.DataFrame({'gene_id':values.index, 'strain':strain, 'mean_expression':values.values}))
    table = richness.merge(pd.concat(frames), on=['gene_id','strain'], validate='one_to_one')
    table['cds_length'] = table.gene_id.map(cds)
    table = table[table.cds_length.gt(0) & table.n_isoforms.le(20) & ~table.gene_id.isin(mt)].copy()
    table['log_expression']=np.log1p(table.mean_expression)
    table['log_cds_length']=np.log(table.cds_length)
    args.output.parent.mkdir(parents=True, exist_ok=True); table.to_csv(args.output,index=False)


if __name__=='__main__': main()
