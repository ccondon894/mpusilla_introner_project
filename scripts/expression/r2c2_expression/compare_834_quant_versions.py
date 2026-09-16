"""Audit 834 quant versions without modifying source files or fitted models."""
import argparse
from pathlib import Path
import hashlib
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--comparison-dir", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--crosswalk", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    root, old, source, out = args.project_root, args.comparison_dir, args.source, args.outdir
    out.mkdir(parents=True, exist_ok=True)
    a = pd.read_csv(source, sep='\t').set_index('Isoform')
    cols = [c for c in a if c.startswith('834_')]
    candidates = sorted(old.glob('*834*clean.quant*')) + sorted(old.glob('new_834*/**/Isoforms.filtered.clean.quant'))
    summary = []
    for path in [source, *candidates]:
        b = pd.read_csv(path, sep='\t').set_index('Isoform')
        assert a.index.is_unique and b.index.is_unique
        ix = a.index.intersection(b.index)
        summary.append(dict(path=str(path), isoforms=len(b), shared=len(ix),
            missing_from_other=len(a.index.difference(b.index)), added_in_other=len(b.index.difference(a.index)),
            gene_changes=int(a.loc[ix, 'Gene'].ne(b.loc[ix, 'Gene']).sum()),
            changed_count_cells=int(a.loc[ix, cols].ne(b.loc[ix, cols]).sum().sum()),
            unique_gene_labels=b.Gene.nunique(), total_counts=int(b[cols].sum().sum()),
            sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    pd.DataFrame(summary).to_csv(out / 'file_comparison.tsv', sep='\t', index=False)
    corrected = pd.read_csv(old / '09092025_834_Isoforms.filtered.clean.quant', sep='\t').set_index('Isoform')
    sq = pd.read_csv(root / 'data/sqanti3_output_834/834_isoforms_classification.filtered.txt', sep='\t').set_index('isoform')
    cross = pd.read_csv(args.crosswalk, sep='\t')
    cross = cross[cross.strain.eq('CCMP1545')].set_index('Isoform')
    diff = pd.DataFrame({'current_quant_gene': a.Gene, 'corrected_quant_gene': corrected.Gene,
                         'model_associated_gene': sq.associated_gene.str.strip(),
                         'model_retained': cross.retained})
    diff['gene_label_changed'] = diff.current_quant_gene.ne(diff.corrected_quant_gene)
    diff['total_reads'] = a[cols].sum(axis=1)
    diff.to_csv(out / 'all_isoform_assignments.tsv', sep='\t', index_label='isoform')
    diff[diff.gene_label_changed].to_csv(out / 'changed_gene_labels.tsv', sep='\t', index_label='isoform')
    sq_rows = []
    for path in sorted((old / 'sqanti3_output_834').glob('834_isoforms_classification.filtered.txt*')):
        b = pd.read_csv(path, sep='\t').set_index('isoform')
        assert set(b.index) == set(sq.index)
        sq_rows.append(dict(path=str(path), **{c + '_differences': int(
            sq[c].fillna('').str.strip().ne(b.loc[sq.index,c].fillna('').str.strip()).sum())
            for c in ['Gene', 'associated_gene']}))
    pd.DataFrame(sq_rows).to_csv(out / 'sqanti_comparison.tsv', sep='\t', index=False)
    changes = diff[diff.gene_label_changed]
    lines = [f'Compared {len(summary)} quant files including the current source.',
        f'Corrected-copy gene-label changes: {len(changes)} / {len(diff)}.',
        f'Changed labels originally starting with Locus: {changes.current_quant_gene.str.startswith("Locus").sum()}.',
        f'Changed-label isoforms currently retained in model: {changes.model_retained.sum()}.',
        'Gene assignments used by the model are reported in all_isoform_assignments.tsv.']
    assert all(x['changed_count_cells'] == 0 and x['shared'] == len(a) and x['isoforms'] == len(a) for x in summary)
    (out / 'summary.txt').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    print(pd.DataFrame(summary)[['path','gene_changes','unique_gene_labels']].to_string(index=False))


if __name__ == '__main__':
    main()
