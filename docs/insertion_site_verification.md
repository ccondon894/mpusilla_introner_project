# Insertion Site Verification Pipeline

This document describes the codon-level introner verification pipeline added on
the `insertion-site-verification` branch. It refines the flank-based ortholog
detection by adding direct evidence from codon/intron position fingerprints
and introner body sequence comparison.

## Motivation

The original ortholog detection (`introner_genotype_matrix_builder.py`) groups
introners across samples by aligning their 100bp flanking sequences. This
approach has two limitations:

1. **No direct evidence for "shared" loci.** Two introners with similar flanks
   are grouped together, but this is no guarantee they are at the same
   biological position. The flank alignment can succeed even when:
   - The introner sequences themselves are completely different (independent
     insertions at the same approximate locus)
   - The introners are at different positions within a gene that has internal
     repeats or tandem-duplicated regions

2. **No way to distinguish ancestral sharing from convergent insertion.** A
   group of introners that all align well across samples could represent:
   - One ancestral insertion event before the common ancestor of all samples
   - Multiple independent insertion events at the same hotspot
   - Annotation/alignment artifacts

The verification pipeline addresses both by computing exact insertion site
fingerprints from splice site annotations + GTF CDS structure, then comparing
these fingerprints across samples within each ortholog group.

## Pipeline overview

```
fix_orientations
    │
    ▼
genotype_matrix.oriented.tsv
    │
    ├──► compute_insertion_fingerprints (Phase 0, per-sample)
    │       │
    │       └──► insertion_fingerprints/{sample}.fingerprints.tsv
    │
    ├──► detect_tandem_duplicates (Phase 0, per-sample)
    │       │
    │       └──► insertion_fingerprints/{sample}.tandems.tsv
    │
    ▼
classify_sharing_status
    │
    ▼
genotype_matrix.verified.tsv      ← adds: sharing_status
    │
    ▼
split_overmerged_orthologs        ← uses fingerprints + tandems + BEDs
    │
    ▼
genotype_matrix.split.tsv         ← adds: original_ortholog_id
    │                                splits over-merged groups
    │                                recovers lost introners from BED
    │
    ▼
reclassify_after_split (rerun classify on split matrix)
    │
    ▼
genotype_matrix.split_verified.tsv
    │
    ▼
compare_introner_sequences        ← extracts bodies, computes pairwise identity
    │
    ▼
genotype_matrix.split_seqverified.tsv
    │                              ← adds: cross_group_identity,
    │                                       within_group_identity,
    │                                       refined_sharing_status
    ▼
annotate_missing_data
    │
    ▼
genotype_matrix.tsv (final)
```

## Scripts (in `scripts/genotyping/`)

### `compute_insertion_fingerprints.py`

For each introner with a gene and splice site annotation, computes the exact
insertion site relative to the host gene's CDS structure.

**Coordinate handling:**
- BED coordinates include 100bp flanks; true body is `[start+100, end-100)`
- For forward-oriented introners: `GT@N` annotation gives the donor splice site
  position directly. Last exonic nucleotide = `body_start + N - 1`.
- For reverse-oriented introners: the GT annotation in `splice_info` is on the
  + strand and is NOT the functional donor. The script reverse-complements the
  body and searches for the GT/GC donor on the - strand, then maps back to
  genomic coordinates.

**Output fields:**
- `location_type`: `cds` if the insertion is in a CDS exon, `intron` if in a
  gap between CDS exons
- `cds_position`, `codon_number`, `codon_offset`, `exon_number`: for `cds`
  insertions, the position relative to the CDS start
- `intron_number`: for `intron` insertions, the 1-indexed intron number in
  transcript order
- `confidence`: `high` if mappable, otherwise an error code (`no_gene`,
  `not_in_cds`, `no_splice_site`, etc.)

**Why intron numbers matter:** miniprot annotates exon boundaries differently
across samples (depending on alignment quality and small indels). This can
shift the codon position of an introner that's actually at the same biological
locus by several codons. The intron number is invariant to these annotation
differences because it's just "which intron of the gene" — independent of the
exact exon lengths.

### `detect_tandem_duplicates.py`

Per-sample preprocessing that finds pairs of introners within a configurable
genomic window (default 10kb) with very high body sequence identity (default
≥95%). Groups them via union-find. Output annotates each introner with
`tandem_cluster_id` and `cluster_size`.

Tandem clusters arise from gene duplications and produce introners that the
flank-based ortholog detection cannot distinguish (because their flanking
regions are also duplicated). The output is used to flag affected ortholog
groups in the splitting step.

### `classify_sharing_status.py`

For each ortholog group, compares fingerprints across presence=1 members using
**chained union-find clustering**. Two members are in the same cluster if
their locus keys are compatible (within tolerance), with transitive merging.
This handles both:
- Continuous noisy distributions (e.g., codons 422-433 with annotation noise)
  → one cluster
- Discrete distinct positions (e.g., codons 282 and 313) → multiple clusters

**Locus key compatibility rules:**
- Both `intron`: same `intron_number`
- Both `cds`: codon difference ≤ tolerance (default 3); offset is NOT required
  to match (small offsets are annotation noise)
- Mixed `cds` and `intron`: cds in exon E is compatible with intron E (right
  after exon E) or intron E−1 (right before exon E). Handles cases where
  miniprot annotated through the introner in some samples but split the gene
  around it in others.

**Sharing status values:**
- `consistent` — within-group only, all members in one chained cluster
- `within_group_discordant` — within-group only, multiple discrete clusters
- `ancestral` — cross-group, exact codon match + same family
- `same_site_diff_family` — exact codon match, different families
- `ambiguous` — cross-group, compatible but not exact (close codons or
  mixed cds/intron)
- `ambiguous_diff_family` — close position, different families
- `independent` — cross-group, codon positions clearly differ
- `uncertain` — insufficient fingerprint data (intergenic, singleton, etc.)

### `split_overmerged_orthologs.py`

Uses chained clustering to detect ortholog groups where Group 1 OR Group 2
contains members at multiple distinct positions (within-group over-merge),
and splits them into sub-groups. Cross-group differences are NOT split —
those represent meaningful biology (independent insertions at the same
approximate locus).

**Lost introner recovery:** When a sub-group is created, the splitter searches
each sample's BED file for additional introners matching the sub-group's
consensus locus. Recovered introners get added as `presence=1` rows. Samples
with no matching BED entry get `presence=2` (genuinely absent). Each sub-group
ends up with a complete sample set with biologically correct calls.

**Output:** Each sub-group gets a new `ortholog_id` like `ortholog_id_4774_split1`,
with the original ID preserved in a new `original_ortholog_id` column.

### `compare_introner_sequences.py`

Extracts introner body sequences from indexed genome FASTAs and computes
pairwise identity using `Bio.Align.PairwiseAligner`. Adds three columns:

- `cross_group_identity`: median G1↔G2 pairwise identity
- `within_group_identity`: representative within-group identity
- `refined_sharing_status`: combines the codon-level status with sequence
  evidence

**Identity calculation note:** Use `alignment.counts()` directly, NOT
`format()` parsing. The format() output includes labels and splits into 60bp
blocks; parsing it manually produces wildly incorrect identity values.

**Threshold rationale:** Two distinct thresholds with different biological
expectations:
- Within-group (default 0.80): within Group 1, ancestral introners are highly
  conserved (~92% median in practice). Catches paralog/mismapping issues.
- Cross-group (default 0.60): accommodates substantial G1↔G2 divergence.
  Used as deciding factor only for `uncertain` (no-codon) cases — codon
  evidence is primary for everything else.

The codon position is the more direct evidence for the insertion event itself.
Sequence identity is largely a measure of introner family relatedness and
isn't a reliable discriminator for ancestral-vs-independent at the locus level.

## Output columns added to the genotype matrix

| Column | Source step | Description |
|---|---|---|
| `sharing_status` | classify_sharing_status | Codon-level classification |
| `original_ortholog_id` | split_overmerged_orthologs | Original flank-based grouping |
| `cross_group_identity` | compare_introner_sequences | G1↔G2 pairwise identity |
| `within_group_identity` | compare_introner_sequences | Within-group identity |
| `refined_sharing_status` | compare_introner_sequences | Final classification |

## Key results

Run on 13-sample matrix with 6,514 original ortholog groups:

| Metric | Original | After verification |
|---|---|---|
| Total ortholog groups | 6,514 | 7,456 |
| Groups split | — | 766 |
| Sub-groups created | — | 942 |
| Lost introners recovered (presence=1) | — | 1,676 |
| G1 fixed (all 11 present) | 4,366 | 3,914 |
| G1 polymorphic (clean) | 761 | **2,095** |
| Cross-group ancestral | 23 | 42 |
| Cross-group independent | 127 | 54 |

The polymorphic count nearly tripled because many "fixed" groups in the
original matrix were actually over-merged groups containing distinct loci that
all had ≥1 sample present. After splitting, the sub-groups have proper
polymorphic patterns.

## Caveats and limitations

### 1. Missing data in polymorphic groups

136 polymorphic groups in Group 1 have at least one sample with `presence=3`
(missing data — the original detection couldn't determine state, typically
due to flank alignment failure). The 547 mixed singletons have similar
issues. **Policy: drop these from allele frequency analysis.**

### 2. Convergent same-site insertions

The `compare_introner_sequences` step identifies a small number (4) of
"`ancestral_low_identity`" cases — exact codon match across G1/G2 but very
low cross-group sequence identity. These may represent real convergent
insertions at the same codon (different introner families landing at the same
hotspot). **Policy: exclude these from Dxy comparisons** since they're not
shared ancestry.

### 3. High-spread within-group groups not split

There are 235 within-group ortholog groups with codon spread >10 that did NOT
get split because they form one chained cluster (continuous noise rather than
discrete positions). Most are likely fine, but a few may contain hidden
over-merges. They are flagged as `within_group_discordant` if the chain
breaks. Worth revisiting if downstream analyses produce suspicious results
in those groups.

### 4. Sequence identity is family-driven, not locus-driven

Empirically, cross-group sequence identity correlates more with introner
family relatedness than with shared insertion ancestry. A pair of family-2
introners shares ~60-80% identity whether or not they came from the same
ancestral insertion. So sequence identity is informational but is NOT used as
the deciding factor for `ancestral` vs `independent` calls — codon position
is. The exception is `uncertain` (no-codon) groups where sequence identity is
the only available evidence.

### 5. Tandem duplications are minor

Despite expectations, only 7 of 766 split groups had tandem cluster members.
Most over-merging is caused by similar flanking regions in distinct loci
within the same gene (gene-internal repeats, near-tandem positions, paralog
mismapping), not by true tandem duplications.

## Snakemake rules added

In `workflow/rules/05_ortholog_detection.smk`:

- `compute_insertion_fingerprints` (Phase 0, per-sample)
- `detect_tandem_duplicates` (Phase 0, per-sample)
- `classify_sharing_status`
- `split_overmerged_orthologs`
- `reclassify_after_split`
- `compare_introner_sequences`

The existing `annotate_missing_data` rule's input was updated to read from
`genotype_matrix.split_seqverified.tsv` (the final output of the verification
pipeline) instead of the original `genotype_matrix.oriented.tsv`.

## Output file locations

In `results/genotyping/`:

- `genotype_matrix.oriented.tsv` — original (Mar 31, preserved unchanged)
- `genotype_matrix.tsv` — original final (Mar 31, preserved unchanged)
- `genotype_matrix.verified.tsv` — after `classify_sharing_status` (new)
- `genotype_matrix.split.tsv` — after `split_overmerged_orthologs` (new)
- `genotype_matrix.split_verified.tsv` — after re-classification on split (new)
- `genotype_matrix.split_seqverified.tsv` — after `compare_introner_sequences`
  (new) — **this is the matrix to use for downstream analysis**

In `results/genotyping/insertion_fingerprints/`:

- `{sample}.fingerprints.tsv` — per-sample codon/intron fingerprints
- `{sample}.tandems.tsv` — per-sample tandem cluster annotations
- `sharing_summary.tsv` — pre-split classification summary
- `sharing_summary_split.tsv` — post-split classification summary
- `sequence_comparison_summary.tsv` — per-group sequence identity summary
- `split_mapping.tsv` — original ortholog ID → split sub-IDs mapping
