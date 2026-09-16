#!/usr/bin/env python3
"""Exon-anchored canonical intron orthology and read-consensus body Dxy.

Exploratory analysis; all coordinates written here are 0-based, half-open.
Run --stage prepare, consensus, or compare; default runs all three.
"""
from __future__ import annotations

import argparse
import concurrent.futures as futures
import hashlib
import importlib.util
import importlib.metadata
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parents[3]
G1 = 'CCMP1545 RCC114 RCC1614 RCC1698 RCC2482 RCC373 RCC465 RCC629 RCC692 RCC693 RCC833'.split()
G2 = ['RCC1749', 'RCC3052']
SAMPLES = G1 + G2
TMP = Path('/scratch1/chris/tmp')
MT = {'CCMP1545#0#scaffold_2': (49807, 1730591),
      'RCC1749#0#intronerless_contig_28': (24999, 2148000)}


def log(msg):
    print(msg, flush=True)


def overlap(a, b, c, d):
    return a < d and c < b


def assembly(sample):
    return ROOT / f'results/assemblies/{sample}.vg_paths.fa'


def annotation(sample):
    return ROOT / f'results/annotations/{sample}.gtf'


def models(sample):
    """Use transcript IDs (called gene_id by the older intron catalog)."""
    raw = defaultdict(list)
    genes = {}
    with annotation(sample).open() as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            f = line.rstrip().split('\t')
            if len(f) != 9 or f[2] != 'exon':
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]+)"', f[8]))
            tx = attrs['transcript_id']
            raw[tx].append((f[0], int(f[3])-1, int(f[4]), f[6]))
            genes[tx] = attrs.get('gene_id', tx)
    result = {}
    with pysam.FastaFile(str(assembly(sample))) as fa:
        for tx, rows in raw.items():
            rows = sorted(set(rows))
            if len({(r[0], r[3]) for r in rows}) != 1:
                continue  # ambiguous multi-contig/strand model
            contig, _, _, strand = rows[0]
            if any(a[2] >= b[1] for a, b in zip(rows, rows[1:])):
                continue  # overlapping or abutting exons, ambiguous structure
            exons = [(r[1], r[2]) for r in rows]
            if strand == '-':
                exons = exons[::-1]
            seqs = [fa.fetch(contig, a, b).upper() for a, b in exons]
            if strand == '-':
                seqs = [str(Seq(s).reverse_complement()) for s in seqs]
            introns = []
            offset = 0
            for i in range(len(exons)-1):
                offset += len(seqs[i])
                a, b = sorted([exons[i], exons[i+1]])
                introns.append(dict(start=a[1], end=b[0], index=i,
                                    junction=offset, left_len=len(seqs[i]),
                                    right_len=len(seqs[i+1])))
            result[tx] = dict(contig=contig, strand=strand, gene_id=genes[tx],
                              start=rows[0][1], end=rows[-1][2],
                              seq=''.join(seqs), introns=introns)
    return result


def junction_map(a, b, mode='global'):
    """Map only boundaries with adjacent ungapped bases on both sides."""
    aligner = PairwiseAligner(mode=mode, match_score=2, mismatch_score=-1,
                             open_gap_score=-5, extend_gap_score=-1)
    alignment = aligner.align(a, b)[0]
    ix = alignment.indices
    columns = np.full(len(a), -1, dtype=int)
    good = ix[0] >= 0
    columns[ix[0, good]] = np.flatnonzero(good)
    mapping = {}
    for j in range(1, len(a)):
        c, d = columns[j-1], columns[j]
        if d == c+1 and ix[1, c] >= 0 and ix[1, d] == ix[1, c]+1:
            mapping[j] = int(ix[1, d])
    return mapping, ix, columns


def anchor_stats(a, b, ix, columns, junction, left, right):
    result = []
    for start, end in [(junction-left, junction), (junction, junction+right)]:
        cols = ix[:, columns[start]:columns[end-1]+1]
        valid = (cols[0] >= 0) & (cols[1] >= 0)
        pairs = [(a[x], b[y]) for x, y in cols[:, valid].T
                 if a[x] in 'ACGT' and b[y] in 'ACGT']
        n = len(pairs)
        result.append((sum(x == y for x, y in pairs)/max(1, cols.shape[1]),
                       n/max(1, end-start), n))
    return result


def match_sample(s, refs, all_models, introner_intervals, old_fixed):
    matched = defaultdict(dict)
    audit = []
    cache = {}
    for n, (lid, (tx, ri)) in enumerate(refs.items(), 1):
        rm = all_models['CCMP1545'][tx]
        sm = all_models[s].get(tx)
        reason = 'PASS'
        si = None
        anchor_identity = 1.0
        if sm is None:
            reason = 'sample_model_missing_or_ambiguous'
        elif s == 'CCMP1545':
            si = ri
        else:
            if tx not in cache:
                mode = 'global' if s in G2 else 'fogsaa'
                forward = junction_map(rm['seq'], sm['seq'], mode)
                reverse = junction_map(sm['seq'], rm['seq'], mode)[0]
                cache[tx] = forward, reverse
            (jm, ix, columns), reverse = cache[tx]
            q = jm.get(ri['junction'])
            options = [i for i in sm['introns'] if i['junction'] == q]
            if len(options) != 1 or reverse.get(q) != ri['junction']:
                reason = 'no_reciprocal_exact_splice_position'
            else:
                si = options[0]
                left = min(40, ri['left_len'], si['left_len'])
                right = min(40, ri['right_len'], si['right_len'])
                stats = anchor_stats(rm['seq'], sm['seq'], ix, columns, ri['junction'], left, right)
                anchor_identity = min(v[0] for v in stats)
                if any(identity < .65 or cov < .8 or n < 20 for identity, cov, n in stats):
                    reason = 'weak_exon_anchor'
        if reason == 'PASS':
            if si['end']-si['start'] < 10:
                reason = 'body_shorter_than_10bp'
            elif any(overlap(si['start'], si['end'], a, b)
                     for a, b in introner_intervals[(s, sm['contig'])]):
                reason = 'current_introner_overlap'
            elif sm['contig'] in MT and overlap(si['start'], si['end'], *MT[sm['contig']]):
                reason = 'direct_MT_overlap'
        audit.append(dict(locus_id=lid, sample=s, status=reason))
        if reason == 'PASS':
            matched[lid][s] = dict(locus_id=lid, sample=s, transcript_id=tx,
                gene_id=rm['gene_id'], contig=sm['contig'], start=si['start'], end=si['end'],
                strand=sm['strand'], intron_index=si['index'], min_anchor_identity=anchor_identity,
                old_matrix_fixed_all13=(tx, ri['index']) in old_fixed)
        if n % 1000 == 0:
            log(f'Orthology {s}: checked {n}/{len(refs)} reference loci')
    return s, dict(matched), audit


def cached_match(sample, refs, all_models, intervals, old_fixed, cache, signature):
    path = cache/f'{sample}.json'
    if path.exists():
        old = json.loads(path.read_text())
        if old['signature'] == signature:
            return old['result']
    result = match_sample(sample, refs, all_models, intervals, old_fixed)
    path.write_text(json.dumps(dict(signature=signature, result=result)))
    return result


def prepare(out, workers):
    catalog = pd.read_csv(ROOT/'results/evolution/non_introner_introns/reference_intron_catalog.tsv', sep='\t')
    matrix = pd.read_csv(ROOT/'results/genotyping/genotype_matrix.final.tsv', sep='\t')
    old = pd.read_csv(ROOT/'results/evolution/non_introner_introns/intron_genotype_matrix.tsv', sep='\t')
    old_fixed = {(r.gene_id, int(r.intron_index)) for r in old.itertuples()
                 if all(getattr(r, s) == 1 for s in SAMPLES)}
    all_models = {}
    for s in SAMPLES:
        all_models[s] = models(s)
    mt_tx = set()
    for s in ['CCMP1545', 'RCC1749']:
        for tx, m in all_models[s].items():
            if m['contig'] in MT and overlap(m['start'], m['end'], *MT[m['contig']]):
                mt_tx.add(tx)
    mt_intervals = []
    for s, mm in all_models.items():
        for tx in sorted(mt_tx & mm.keys()):
            m = mm[tx]
            mt_intervals.append(dict(sample=s, transcript_id=tx, contig=m['contig'],
                                     start=m['start'], end=m['end']))
    pd.DataFrame(mt_intervals).to_csv(out/'mt_orthologous_gene_intervals.tsv', sep='\t', index=False)
    introner_intervals = defaultdict(list)
    for r in matrix[matrix.presence == 1].itertuples():
        introner_intervals[(r.sample, r.contig)].append((int(r.start)+100, int(r.end)-100))

    refs = {}
    audit = []
    for r in catalog.itertuples():
        tx = r.gene_id
        lid = f'{tx}__intron_{r.intron_idx}'
        m = all_models['CCMP1545'].get(tx)
        reason = None
        if tx in mt_tx:
            reason = 'MT_gene_in_either_reference'
        elif m is None:
            reason = 'reference_model_missing_or_ambiguous'
        else:
            found = [i for i in m['introns'] if i['start'] == r.start-1 and i['end'] == r.end]
            if len(found) != 1:
                reason = 'catalog_annotation_mismatch'
            elif any(overlap(found[0]['start'], found[0]['end'], a, b)
                     for a, b in introner_intervals[('CCMP1545', m['contig'])]):
                reason = 'current_introner_overlap'
        if reason:
            audit.append(dict(locus_id=lid, sample='CCMP1545', status=reason))
        else:
            refs[lid] = (tx, found[0])
    # No intron-body sequence similarity or old presence call is used for entry.
    matched = defaultdict(dict)
    cache = out/'orthology_cache'
    cache.mkdir(exist_ok=True)
    sources = [annotation(s) for s in SAMPLES] + [assembly(s) for s in SAMPLES]
    sources += [ROOT/'results/genotyping/genotype_matrix.final.tsv',
                ROOT/'results/evolution/non_introner_introns/reference_intron_catalog.tsv',
                ROOT/'results/evolution/non_introner_introns/intron_genotype_matrix.tsv']
    signature = hashlib.sha256(json.dumps(dict(method='exon_orthology_v2_global_G2',
                   inputs=[fingerprint(p) for p in sources], mt=MT),sort_keys=True).encode()).hexdigest()
    with futures.ProcessPoolExecutor(max_workers=workers) as pool:
        jobs = [pool.submit(cached_match, s, refs, all_models, introner_intervals, old_fixed, cache, signature) for s in G2+G1]
        for job in futures.as_completed(jobs):
            s, hits, sample_audit = job.result()
            audit.extend(sample_audit)
            for lid, values in hits.items():
                matched[lid].update(values)
            log(f'Orthology {s}: {len(hits)} matched introns')
    rows = [v[s] for v in matched.values() if len(v) == len(SAMPLES) for s in SAMPLES]
    df = pd.DataFrame(rows)
    if df.empty:
        raise RuntimeError('No complete orthologous canonical intron loci')
    # Deduplicate identical genomic reference intervals across transcript isoforms.
    refrows = df[df['sample'] == 'CCMP1545'].sort_values('locus_id')
    keep = set(refrows.drop_duplicates(['contig', 'start', 'end']).locus_id)
    df = df[df.locus_id.isin(keep)]
    # Exclude loci mapping to the same physical intron in any other sample.
    ambiguous = df[df.duplicated(['sample', 'contig', 'start', 'end'], keep=False)].locus_id.unique()
    df = df[~df.locus_id.isin(ambiguous)]
    df.to_csv(out/'orthologous_introns.tsv', sep='\t', index=False)
    pd.DataFrame(audit).to_csv(out/'orthology_audit.tsv', sep='\t', index=False)
    for s in SAMPLES:
        sub = df[df['sample'] == s].sort_values(['contig', 'start', 'end'])
        with (out/f'{s}.bed').open('w') as h:
            for r in sub.itertuples():
                h.write(f'{r.contig}\t{r.start}\t{r.end}\t{r.locus_id}\t0\t{r.strand}\n')
    summary = dict(reference_catalog_loci=len(catalog), old_matrix_fixed_all13=len(old_fixed),
                   MT_transcripts=len(mt_tx), complete_unique_exon_anchored_loci=df.locus_id.nunique(),
                   genes=df.gene_id.nunique(), old_matrix_fixed_loci_retained=int(df[df['sample']=='CCMP1545'].old_matrix_fixed_all13.sum()),
                   duplicate_target_loci_removed=len(ambiguous))
    (out/'preparation_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    log(summary)


def run(cmd, logfile, stdout=None):
    logfile.write(' '.join(map(str, cmd))+'\n')
    logfile.flush()
    subprocess.run(list(map(str, cmd)), check=True,
                   stdout=stdout if stdout is not None else logfile, stderr=logfile)


def fingerprint(path, content=False):
    path = Path(path)
    st = path.stat()
    d = dict(path=str(path.resolve()), size=st.st_size, mtime_ns=st.st_mtime_ns)
    if content:
        d['sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    return d


def consensus_sample(sample, out):
    """Reuse workflow caller and QC; extract regional consensus before indels shift coordinates."""
    bed = out/f'{sample}.bed'
    fa = assembly(sample)
    bam = ROOT/f'results/genotyping/coverage/bams/{sample}.sorted.bam'
    qc_script = ROOT/'scripts/evolution/qc_introner_body_consensus.py'
    old_qc = ROOT/f'results/evolution/introner_body_consensus/{sample}.introner_body.qc.tsv'
    provenance = dict(method='canonical_regional_consensus_v1',
                      bed_sha256=hashlib.sha256(bed.read_bytes()).hexdigest(),
                      assembly=fingerprint(fa), qc_script=fingerprint(qc_script, True),
                      original_qc=fingerprint(old_qc, True))
    if sample != 'CCMP1545':
        provenance['bam'] = fingerprint(bam)
    done = out/f'{sample}.complete.json'
    cons = out/f'{sample}.consensus.fa'
    qc_path = out/f'{sample}.qc.tsv'
    if done.exists() and json.loads(done.read_text()) == provenance and cons.exists() and qc_path.exists():
        log(f'Consensus {sample}: using completed result')
        return
    regions = []
    with bed.open() as h:
        for line in h:
            chrom, start, end, name, _, strand = line.rstrip().split('\t')
            regions.append((chrom, int(start), int(end), name, strand))
    with tempfile.TemporaryDirectory(prefix=f'canonical_dxy_{sample}_', dir=TMP) as tmpname, \
            (out/f'{sample}.consensus.log').open('w') as logfile:
        tmp = Path(tmpname)
        ref_regions = tmp/'regions.fa'
        with pysam.FastaFile(str(fa)) as fh, ref_regions.open('w') as h:
            for chrom, start, end, name, strand in regions:
                h.write(f'>{chrom}:{start+1}-{end}\n{fh.fetch(chrom, start, end)}\n')
        if sample == 'CCMP1545':
            result_regions = ref_regions
            qb = []
        else:
            log(f'Consensus {sample}: fetching reads at {len(regions)} introns')
            subset = tmp/'reads.bam'
            run(['samtools', 'view', '-M', '-L', bed, '-b', '-o', subset, bam], logfile)
            run(['samtools', 'index', subset], logfile)
            pileup = tmp/'pileup.bcf'
            run(['bcftools', 'mpileup', '-d', '100', '-q', '30', '-Q', '20', '-A',
                 '-f', fa, '-R', bed, '-Ob', '-o', pileup, subset], logfile)
            calls = tmp/'calls.bcf'
            run(['bcftools', 'call', '-m', '--ploidy', '1', '-Ob', '-o', calls, pileup], logfile)
            norm = tmp/'norm.vcf.gz'
            run(['bcftools', 'norm', '-f', fa, '-m', '+both', '-Oz', '-o', norm, calls], logfile)
            filt = out/f'{sample}.filtered.vcf.gz'
            run(['bcftools', 'filter', '-i', 'DP>=10', '-Oz', '-o', filt, norm], logfile)
            run(['bcftools', 'index', '-f', '-t', filt], logfile)
            result_regions = tmp/'consensus_regions.fa'
            run(['bcftools', 'consensus', '-a', 'N', '-f', ref_regions, '-o', result_regions, filt], logfile)
            qb = ['--bam', subset]
        seqs = {r.id: str(r.seq).upper() for r in SeqIO.parse(result_regions, 'fasta')}
        with cons.open('w') as h:
            for chrom, start, end, name, strand in regions:
                key = f'{chrom}:{start+1}-{end}'
                h.write(f'>{name}\n{seqs[key]}\n')
        log(f'Consensus {sample}: applying existing introner-body QC')
        run([sys.executable, qc_script, '--sample', sample, '--bed', bed,
             '--consensus-fasta', cons, '--output', qc_path,
             '--reference-fasta', fa, *qb], logfile)
        # Use the SAME per-sample depth baseline as the existing introner QC.
        # Preserve the canonical-only ratio as an audit column.
        qc = pd.read_csv(qc_path, sep='\t')
        original = pd.read_csv(old_qc, sep='\t')
        baseline = float(original.sample_median_body_depth.iloc[0])
        qc['canonical_median_body_depth'] = qc.sample_median_body_depth
        qc['canonical_depth_ratio'] = qc.depth_ratio
        if sample != 'CCMP1545':
            if not np.isfinite(baseline) or baseline <= 0:
                raise ValueError(f'Invalid original depth baseline for {sample}')
            qc['sample_median_body_depth'] = baseline
            qc['depth_ratio'] = qc.mean_depth / baseline
        spec = importlib.util.spec_from_file_location('body_qc', qc_script)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        qc['qc_fail_reason'] = [module.fail_reasons(r) for r in qc.to_dict('records')]
        qc['qc_pass'] = qc.qc_fail_reason == 'PASS'
        qc.to_csv(qc_path, sep='\t', index=False)
        log(f'Consensus {sample}: {int(qc.qc_pass.sum())}/{len(qc)} pass')
    done.write_text(json.dumps(provenance, indent=2)+'\n')


def consensus(out, workers):
    with futures.ThreadPoolExecutor(max_workers=workers) as pool:
        jobs = [pool.submit(consensus_sample, s, out) for s in SAMPLES]
        for job in futures.as_completed(jobs):
            job.result()
    qc = pd.concat([pd.read_csv(out/f'{s}.qc.tsv', sep='\t') for s in SAMPLES], ignore_index=True)
    qc.to_csv(out/'canonical_consensus_qc.tsv', sep='\t', index=False)


def dxy_counts(seqs, trim=0):
    """Ratio of summed differences to summed callable sites over all 22 pairs."""
    arrays = {}
    for s in SAMPLES:
        a = np.frombuffer(seqs[s].upper().encode(), dtype='S1').copy()
        if trim:
            bases = np.flatnonzero(a != b'-')
            a[bases[:trim]] = b'N'
            a[bases[-trim:]] = b'N'
        arrays[s] = a
    if len({len(a) for a in arrays.values()}) != 1:
        raise ValueError('Unequal alignment lengths')
    pairs = []
    for s in G1:
        for t in G2:
            a, b = arrays[s], arrays[t]
            valid = np.isin(a, [b'A', b'C', b'G', b'T']) & np.isin(b, [b'A', b'C', b'G', b'T'])
            n = int(valid.sum())
            diff = int(((a != b) & valid).sum())
            shorter = min(int(np.isin(a, [b'A', b'C', b'G', b'T']).sum()),
                          int(np.isin(b, [b'A', b'C', b'G', b'T']).sum()))
            pairs.append(dict(sample_group1=s, sample_group2=t, differences=diff,
                              callable_pair_sites=n, aligned_fraction_of_shorter=n/max(1, shorter)))
    sites = sum(r['callable_pair_sites'] for r in pairs)
    diffs = sum(r['differences'] for r in pairs)
    return dict(dxy=diffs/sites if sites else np.nan, differences=diffs,
                callable_pair_sites=sites, min_pair_sites=min(r['callable_pair_sites'] for r in pairs),
                min_aligned_fraction=min(r['aligned_fraction_of_shorter'] for r in pairs)), pairs


def align_canonical(lid, sequences, out):
    dest = out/'alignments'/f'{lid}.mafft.fa'
    digest = hashlib.sha256(json.dumps(sequences, sort_keys=True).encode()).hexdigest()
    stamp = dest.with_suffix('.sha256')
    if not (dest.exists() and stamp.exists() and stamp.read_text().strip() == digest):
        with tempfile.TemporaryDirectory(prefix='canonical_mafft_', dir=TMP) as td:
            fa = Path(td)/'input.fa'
            fa.write_text(''.join(f'>{s}\n{sequences[s]}\n' for s in SAMPLES))
            cp = subprocess.run(['mafft', '--thread', '1', '--adjustdirection', '--maxiterate',
                                 '1000', '--globalpair', '--quiet', str(fa)],
                                capture_output=True, text=True, timeout=300, check=True,
                                env={**os.environ, 'TMPDIR': str(TMP)})
            dest.write_text(cp.stdout)
            stamp.write_text(digest+'\n')
    seqs = {r.id.removeprefix('_R_'): str(r.seq) for r in SeqIO.parse(dest, 'fasta')}
    if set(seqs) != set(SAMPLES):
        raise ValueError(f'{lid}: incomplete alignment')
    metrics, pairs = dxy_counts(seqs)
    metrics['dxy_trim2'] = dxy_counts(seqs, trim=2)[0]['dxy']
    return lid, metrics, pairs


def introner_data(out, use_regional=True):
    body = pd.read_csv(ROOT/'results/evolution/diversity_metrics/shared_introner_body_dxy_200bp.tsv', sep='\t')
    body = body[body.ancestry_class == 'ancestral']
    if not use_regional:
        body = body[body.dxy_introner.notna()]
    matrix = pd.read_csv(ROOT/'results/genotyping/genotype_matrix.final.tsv', sep='\t')
    existing_qc = pd.read_csv(ROOT/'results/evolution/introner_body_consensus/introner_body_consensus.qc.tsv', sep='\t').set_index(['sample','ortholog_id'])
    regional_dir = out/'introner_regional_sensitivity'
    if use_regional:
        regional = pd.read_csv(regional_dir/'legacy_vs_regional_dxy.tsv',sep='\t').set_index('ortholog_id')
        rqc = pd.read_csv(regional_dir/'regional_consensus_qc.tsv',sep='\t')
        nr = rqc[rqc.qc_pass.astype(str).str.lower().eq('true')].groupby('ortholog_id')['sample'].nunique()
        regional_keep = set(nr[nr==13].index)
    mt = pd.read_csv(out/'mt_orthologous_gene_intervals.tsv', sep='\t')
    mt_index = defaultdict(list)
    for r in mt.itertuples():
        mt_index[(r.sample, r.contig)].append((r.start, r.end))
    ref_models = models('CCMP1545')
    rows, pairs, excluded = [], [], []
    for r in body.itertuples():
        lid = r.ortholog_id
        loc = matrix[(matrix.ortholog_id == lid) & (matrix.presence == 1) & matrix['sample'].isin(SAMPLES)]
        if len(loc) != 13 or loc['sample'].nunique() != 13:
            excluded.append(dict(locus_id=lid, reason='current_matrix_not_fixed_all13'))
            continue
        if not set(loc.cross_group_status).issubset({'ancestral','likely_ancestral','ancestral_low_identity'}):
            excluded.append(dict(locus_id=lid,reason='current_origin_not_ancestral'))
            continue
        for x in loc.itertuples():
            q = existing_qc.loc[(x.sample, lid)]
            if (q.contig != x.contig or q.body_start != x.start+100 or q.body_end != x.end-100
                    or (not use_regional and str(q.qc_pass).lower() != 'true')):
                raise ValueError(f'{lid}/{x.sample}: existing introner QC is stale or failed')
        if any((x.contig in MT and overlap(x.start+100, x.end-100, *MT[x.contig])) or
               any(overlap(x.start+100, x.end-100, a, b) for a, b in mt_index[(x.sample, x.contig)])
               for x in loc.itertuples()):
            excluded.append(dict(locus_id=lid, reason='MT_or_orthologous_MT_gene'))
            continue
        if use_regional and lid not in regional_keep:
            excluded.append(dict(locus_id=lid, reason='regional_consensus_failed_QC'))
            continue
        if pd.notna(r.dxy_introner):
            fa = ROOT/f'results/evolution/alignments/introner_body/{lid}.introner_body.mafft.fa'
            seqs = {x.id.removeprefix('_R_'): str(x.seq) for x in SeqIO.parse(fa, 'fasta')}
            if set(seqs) != set(SAMPLES):
                raise ValueError(f'{lid}: existing alignment missing samples')
            metrics, pp = dxy_counts(seqs)
            if not np.isclose(metrics['dxy'], r.dxy_introner, atol=1e-12):
                raise ValueError(f'{lid}: recomputed Dxy does not reproduce existing table')
        if use_regional:
            fa = regional_dir/'alignments'/f'{lid}.mafft.fa'
            seqs = {x.id.removeprefix('_R_'): str(x.seq) for x in SeqIO.parse(fa, 'fasta')}
            metrics, pp = dxy_counts(seqs)
            if not np.isclose(metrics['dxy'],regional.loc[lid,'dxy'],atol=1e-12):
                raise ValueError(f'{lid}: regional Dxy does not reproduce sensitivity table')
        ref = loc[loc['sample']=='CCMP1545'].iloc[0]
        start, end = int(ref.start)+100, int(ref.end)-100
        genes = {m['gene_id'] for m in ref_models.values() if m['contig'] == ref.contig
                 and overlap(start, end, m['start'], m['end'])}
        gene = next(iter(genes)) if len(genes) == 1 else f'unassigned:{lid}'
        lengths = [len(seqs[s].replace('-', '')) for s in SAMPLES]
        rows.append(dict(locus_id=lid, region_class='Introner', gene_id=gene,
                         gene_assignment='unique_overlap' if len(genes)==1 else 'unassigned_or_ambiguous',
                         ref_body_length=end-start, median_body_length=float(np.median(lengths)),
                         legacy_dxy=r.dxy_introner, dxy_trim2=dxy_counts(seqs, trim=2)[0]['dxy'], **metrics))
        pairs.extend(dict(locus_id=lid, region_class='Introner', **x) for x in pp)
    name = 'introner_exclusions.tsv' if use_regional else 'legacy_introner_exclusions.tsv'
    pd.DataFrame(excluded, columns=['locus_id','reason']).to_csv(out/name, sep='\t', index=False)
    return rows, pairs


def compare(out, workers):
    subprocess.run([sys.executable, Path(__file__).with_name('check_introner_regional_consensus.py'),
                    '--output',out/'introner_regional_sensitivity'],check=True)
    orth = pd.read_csv(out/'orthologous_introns.tsv', sep='\t')
    qc = pd.read_csv(out/'canonical_consensus_qc.tsv', sep='\t')
    passed = qc[qc.qc_pass.astype(str).str.lower().eq('true')].groupby('ortholog_id')['sample'].nunique()
    good = set(passed[passed == 13].index)
    ref = orth[orth['sample']=='CCMP1545'].set_index('locus_id')
    seqs = {}
    for s in SAMPLES:
        strands = orth[orth['sample']==s].set_index('locus_id').strand.to_dict()
        for record in SeqIO.parse(out/f'{s}.consensus.fa', 'fasta'):
            if record.id not in good:
                continue
            seq = str(record.seq)
            if strands[record.id] == '-':
                seq = str(Seq(seq).reverse_complement())
            seqs.setdefault(record.id, {})[s] = seq
    (out/'alignments').mkdir(exist_ok=True)
    rows, pairs = introner_data(out)
    legacy_rows, _ = introner_data(out,use_regional=False)
    pd.DataFrame(legacy_rows).to_csv(out/'legacy_introner_body_dxy.tsv',sep='\t',index=False)
    log(f'Comparing {len(rows)} ancestral introners with {len(good)} QC-passing canonical introns')
    with futures.ThreadPoolExecutor(max_workers=workers) as pool:
        jobs = [pool.submit(align_canonical, lid, sequences, out) for lid, sequences in sorted(seqs.items())]
        for n, job in enumerate(futures.as_completed(jobs), 1):
            lid, metrics, pp = job.result()
            rr = ref.loc[lid]
            lengths = [len(s) for s in seqs[lid].values()]
            rows.append(dict(locus_id=lid, region_class='Canonical intron', gene_id=rr.gene_id,
                             gene_assignment='annotated_transcript', ref_body_length=int(rr.end-rr.start),
                             median_body_length=float(np.median(lengths)),
                             old_matrix_fixed_all13=rr.old_matrix_fixed_all13, **metrics))
            pairs.extend(dict(locus_id=lid, region_class='Canonical intron', **x) for x in pp)
            if n % 100 == 0 or n == len(jobs):
                log(f'MAFFT and Dxy: {n}/{len(jobs)} canonical loci')
    df = pd.DataFrame(rows).sort_values(['region_class','locus_id'])
    df['alignment_qc_pass'] = df.min_pair_sites > 0
    df.to_csv(out/'body_dxy_per_locus.tsv', sep='\t', index=False)
    pd.DataFrame(pairs).to_csv(out/'body_dxy_per_sample_pair.tsv', sep='\t', index=False)
    summarize(df[df.alignment_qc_pass].copy(), out)


def summarize(df, out):
    from scipy.stats import mannwhitneyu, wilcoxon
    from scipy.optimize import linear_sum_assignment
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(20260910)
    it = df[df.region_class=='Introner'].copy()
    ca = df[df.region_class=='Canonical intron'].copy()
    if it.empty or ca.empty:
        raise ValueError('Both body classes are required')
    tests = []
    def mwu(a, b, name):
        if not len(a) or not len(b):
            return
        stat = mannwhitneyu(a, b, alternative='two-sided')
        tests.append(dict(comparison=name, test='Mann-Whitney U (two-sided)', n_introner=len(a),
                          n_control=len(b), mean_introner=float(np.mean(a)), mean_control=float(np.mean(b)),
                          mean_difference=float(np.mean(a)-np.mean(b)), statistic=float(stat.statistic), pvalue=float(stat.pvalue)))
    mwu(it.dxy, ca.dxy, 'all_callable_bodies')
    legacy = pd.read_csv(out/'legacy_introner_body_dxy.tsv',sep='\t')
    mwu(legacy.dxy,ca.dxy,'legacy_genome_wide_introner_consensus')
    mwu(it.dxy_trim2.dropna(), ca.dxy_trim2.dropna(), 'remove_2bp_at_each_splice_end')
    hi = df[(df.min_pair_sites >= 50) & (df.min_aligned_fraction >= .5)]
    mwu(hi[hi.region_class=='Introner'].dxy, hi[hi.region_class=='Canonical intron'].dxy,
        'minimum_50_sites_and_half_shorter_sequence_each_pair')
    ca_old = ca[ca.old_matrix_fixed_all13.astype(str).str.lower().eq('true')]
    mwu(it.dxy, ca_old.dxy, 'old_similarity_based_presence_calls_subset')
    # Match without replacement on reference body length; report 2-fold caliper failures.
    cost = np.abs(np.log(it.ref_body_length.to_numpy()[:,None] / ca.ref_body_length.to_numpy()[None,:]))
    ii, jj = linear_sum_assignment(cost)
    matches = []
    for i,j in zip(ii,jj):
        a, b = it.iloc[i], ca.iloc[j]
        matches.append(dict(introner_locus=a.locus_id, canonical_locus=b.locus_id,
                            introner_gene=a.gene_id, canonical_gene=b.gene_id,
                            introner_length=a.ref_body_length, canonical_length=b.ref_body_length,
                            length_ratio=a.ref_body_length/b.ref_body_length,
                            within_2fold=bool(cost[i,j] <= np.log(2)),
                            introner_dxy=a.dxy, canonical_dxy=b.dxy, difference=a.dxy-b.dxy))
    match = pd.DataFrame(matches)
    match.to_csv(out/'length_matched_pairs.tsv', sep='\t', index=False)
    match = match[match.within_2fold]
    def paired(a, b, name):
        if not len(a):
            return
        diff = np.asarray(a)-np.asarray(b)
        stat = wilcoxon(diff, alternative='two-sided') if np.any(diff) else None
        tests.append(dict(comparison=name, test='Wilcoxon signed-rank (two-sided)',
                          n_introner=len(a), n_control=len(b), mean_introner=float(np.mean(a)),
                          mean_control=float(np.mean(b)), mean_difference=float(np.mean(diff)),
                          statistic=float(stat.statistic) if stat else 0., pvalue=float(stat.pvalue) if stat else 1.))
    paired(match.introner_dxy, match.canonical_dxy, 'length_matched_without_replacement_2fold_caliper')
    # Each shared gene contributes one introner mean and one canonical mean.
    gm = df[~df.gene_id.str.startswith('unassigned:')].groupby(['gene_id','region_class']).dxy.mean().unstack().dropna()
    gm['difference'] = gm['Introner']-gm['Canonical intron']
    gm.to_csv(out/'same_gene_comparison.tsv', sep='\t')
    paired(gm['Introner'], gm['Canonical intron'], 'same_gene_mean_bodies')
    # Joint gene-block bootstrap retains dependence within and between body classes.
    grouped = df.groupby(['gene_id','region_class']).dxy.agg(['sum','count']).unstack(fill_value=0)
    sums_i = grouped[('sum','Introner')].to_numpy()
    sums_c = grouped[('sum','Canonical intron')].to_numpy()
    ns_i = grouped[('count','Introner')].to_numpy()
    ns_c = grouped[('count','Canonical intron')].to_numpy()
    boot = []
    for _ in range(10000):
        idx = rng.integers(0,len(grouped),len(grouped))
        if ns_i[idx].sum() and ns_c[idx].sum():
            boot.append(sums_i[idx].sum()/ns_i[idx].sum()-sums_c[idx].sum()/ns_c[idx].sum())
    ci = np.quantile(boot,[.025,.975])
    pd.DataFrame(tests).to_csv(out/'statistical_tests.tsv', sep='\t', index=False)
    summaries = df.groupby('region_class').agg(n=('dxy','size'), gene_blocks=('gene_id','nunique'),
                   mean_dxy=('dxy','mean'), median_dxy=('dxy','median'), sd_dxy=('dxy','std'),
                   mean_ref_length=('ref_body_length','mean'), median_ref_length=('ref_body_length','median'))
    assigned = ~df.gene_id.str.startswith('unassigned:')
    summaries['annotated_genes'] = df[assigned].groupby('region_class').gene_id.nunique().reindex(summaries.index,fill_value=0)
    summaries['unassigned_loci'] = df[~assigned].groupby('region_class').size().reindex(summaries.index,fill_value=0)
    summaries.to_csv(out/'class_summary.tsv', sep='\t')
    # Recheck body/flank contrast on the same retained introner loci.
    flanks = pd.read_csv(ROOT/'results/evolution/diversity_metrics/all_samples_diversity_metrics_200bp.tsv', sep='\t')
    paired_flanks = it.merge(flanks[['ortholog_id','dxy_group1_group2']], left_on='locus_id', right_on='ortholog_id').dropna(subset=['dxy_group1_group2'])
    paired_flanks.to_csv(out/'introner_body_vs_paired_flanks.tsv', sep='\t', index=False)
    paired(paired_flanks.dxy, paired_flanks.dxy_group1_group2, 'introner_body_vs_same_locus_flanks')
    pd.DataFrame(tests).to_csv(out/'statistical_tests.tsv', sep='\t', index=False)
    fig, axes = plt.subplots(1,3,figsize=(12,4.2), layout='constrained')
    colors = ['#bc272d','#0b81a2']
    groups = [ca.dxy.to_numpy(),it.dxy.to_numpy()]
    vp = axes[0].violinplot(groups,showmeans=True,showextrema=False)
    for artist,col in zip(vp['bodies'],colors):
        artist.set_facecolor(col); artist.set_alpha(.65)
    axes[0].boxplot(groups, widths=.12,showfliers=False)
    axes[0].set_xticks([1,2],[f'Canonical introns\nn = {len(ca):,}',f'Ancestral introners\nn = {len(it):,}'])
    axes[0].set_ylabel('Body Dxy'); axes[0].set_title('All callable bodies')
    for _,r in match.iterrows():
        axes[1].plot([1,2],[r.canonical_dxy,r.introner_dxy],color='#999999',alpha=.45,lw=.8)
    axes[1].scatter(np.ones(len(match)),match.canonical_dxy,color=colors[0],s=14)
    axes[1].scatter(np.full(len(match),2),match.introner_dxy,color=colors[1],s=14)
    axes[1].set_xticks([1,2],['Canonical introns','Introners']); axes[1].set_title(f'Length matched ({len(match)} pairs)')
    for _,r in gm.iterrows():
        axes[2].plot([1,2],[r['Canonical intron'],r['Introner']],color='#999999',alpha=.6,lw=.8)
    axes[2].scatter(np.ones(len(gm)),gm['Canonical intron'],color=colors[0],s=18)
    axes[2].scatter(np.full(len(gm),2),gm['Introner'],color=colors[1],s=18)
    axes[2].set_xticks([1,2],['Canonical introns','Introners']); axes[2].set_title(f'Same gene ({len(gm)} genes)')
    for ax in axes:
        ax.spines[['top','right']].set_visible(False); ax.set_ylim(bottom=0)
    test_lookup = {r['comparison']: r for r in tests}
    for ax, name in zip(axes,['all_callable_bodies','length_matched_without_replacement_2fold_caliper','same_gene_mean_bodies']):
        if name in test_lookup:
            ax.set_title(ax.get_title()+f"\np = {test_lookup[name]['pvalue']:.3g}")
    fig.savefig(out/'canonical_vs_introner_body_dxy.png',dpi=220)
    plt.close(fig)
    summary = dict(mean_introner=float(it.dxy.mean()), mean_canonical=float(ca.dxy.mean()),
                   ratio=float(it.dxy.mean()/ca.dxy.mean()), mean_difference=float(it.dxy.mean()-ca.dxy.mean()),
                   gene_block_bootstrap_95ci_difference=ci.tolist(), bootstrap_replicates=len(boot), seed=20260910,
                   n_introner=len(it), n_canonical=len(ca), same_gene_pairs=len(gm), length_pairs=len(match),
                   same_locus_flank_mean=float(paired_flanks.dxy_group1_group2.mean()))
    (out/'comparison_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    log(json.dumps(summary,indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--stage', choices=['all', 'prepare', 'consensus', 'compare'], default='all')
    p.add_argument('--workers', type=int, default=4)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    TMP.mkdir(parents=True, exist_ok=True)
    os.environ['TMPDIR'] = str(TMP)
    tempfile.tempdir = str(TMP)
    inputs = [ROOT/'results/genotyping/genotype_matrix.final.tsv',
              ROOT/'results/evolution/non_introner_introns/reference_intron_catalog.tsv',
              ROOT/'results/evolution/non_introner_introns/intron_genotype_matrix.tsv',
              ROOT/'results/evolution/diversity_metrics/shared_introner_body_dxy_200bp.tsv',
              ROOT/'results/evolution/diversity_metrics/all_samples_diversity_metrics_200bp.tsv',
              ROOT/'scripts/evolution/qc_introner_body_consensus.py', Path(__file__).resolve(),
              Path(__file__).with_name('check_introner_regional_consensus.py'),
              *[annotation(s) for s in SAMPLES]]
    provenance = dict(command=sys.argv, python=sys.version, inputs=[fingerprint(x, True) for x in inputs],
                      assemblies=[fingerprint(assembly(s)) for s in SAMPLES],
                      groups=dict(group1=G1,group2=G2), MT_intervals_0based=MT,
                      package_versions={x: importlib.metadata.version(x) for x in ['biopython','pysam','numpy','pandas','scipy','matplotlib']})
    provenance['tool_versions'] = {}
    for tool in ['mafft','samtools','bcftools']:
        cp = subprocess.run([tool,'--version'],capture_output=True,text=True,errors='replace',check=True)
        provenance['tool_versions'][tool] = (cp.stdout+cp.stderr).strip().splitlines()[0]
    (args.output/f'provenance_{args.stage}.json').write_text(json.dumps(provenance,indent=2)+'\n')
    if args.stage in ['all', 'prepare']:
        prepare(args.output, args.workers)
    if args.stage in ['all', 'consensus']:
        consensus(args.output, args.workers)
    if args.stage in ['all', 'compare']:
        compare(args.output, args.workers)


if __name__ == '__main__':
    main()
