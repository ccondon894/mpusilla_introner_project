# Sharing Status Categories

This document describes the classification columns in the genotype matrix
(`genotype_matrix.tsv`) and the meaning of each category.

## Overview

Each ortholog group receives two independent classifications:

1. **`within_group_status`** — consistency of insertion site within each
   clade. Use this to filter groups for within-group nucleotide diversity
   (pi) calculations.

2. **`cross_group_status`** — ancestry relationship between G1 and G2.
   Use this to identify ancestral vs independent introners for between-group
   divergence (dxy) calculations.

Both columns are initially set by `classify_sharing_status.py` using
codon-level insertion fingerprints, then refined by
`compare_introner_sequences.py` using introner body sequence identity.
The values in the final matrix reflect the refined classification.

Two audit columns record the sequence evidence:
- **`within_group_identity`**: pairwise identity between two members of
  the largest clade (float, empty if not computed)
- **`cross_group_identity`**: median G1-vs-G2 pairwise identity
  (float, empty if not computed)

### Context: Group 1 vs Group 2

The 13 *M. pusilla* strains fall into two phylogenetic clades:

- **Group 1 (G1)**: CCMP1545, RCC114, RCC1614, RCC1698, RCC2482, RCC373,
  RCC465, RCC629, RCC692, RCC693, RCC833 (11 strains, introner-rich)
- **Group 2 (G2)**: RCC1749, RCC3052 (2 strains, introner-depleted)

---

## `within_group_status`

| Category | Meaning | Use for pi? |
|---|---|---|
| **consistent** | All members within each clade share the same insertion site. Verified by either codon fingerprints (genic loci) or body sequence identity >= 80% with same family (intergenic loci). | Yes |
| **singleton** | Only 1 presence=1 member in the entire group. No comparison possible. | No (n=1) |
| **low_identity** | Body sequence identity < 80% between members of the same clade. May indicate paralog mismerging, rapid introner evolution, or partial deletions. | Use with caution |
| **discordant** | Members within a single clade cluster to multiple distinct insertion sites. Rare after splitting. | No |
| **uncertain** | Insufficient data (could not extract body sequences). | No |

### How within-group classification works

1. **Codon fingerprinting** (`classify_sharing_status.py`): For genic loci,
   compares the hybrid aa-context + codon/intron position within each clade.
   If all members cluster together -> `consistent`. If multiple clusters ->
   `discordant`. If insufficient fingerprints -> `uncertain`.

2. **Sequence refinement** (`compare_introner_sequences.py`): For groups
   still `uncertain` (mostly intergenic), extracts introner body sequences
   and checks: if within-group identity >= 80% AND all members share the
   same introner family -> reclassified as `consistent`. If identity < 80%
   -> reclassified as `low_identity`.

---

## `cross_group_status`

| Category | Meaning | Use for dxy? |
|---|---|---|
| **ancestral** | Same insertion site (exact or aa-context match) + same family. Strong evidence the introner was present before G1/G2 divergence. | Yes (shared) |
| **ancestral_low_identity** | Same as ancestral, but cross-group body sequence identity < 60%. Codon evidence says ancestral but sequences are very divergent. | Review case-by-case |
| **likely_ancestral** | No codon data (intergenic). Cross-group body sequence identity >= 60%. | Yes (shared, lower confidence) |
| **independent** | Different insertion site, different families, or only legacy codon number match without conserved protein context around the insertion site. | Yes (independent) |
| **likely_independent** | No codon data (intergenic). Cross-group body sequence identity < 60%. | Yes (independent, lower confidence) |
| **uncertain** | Insufficient data for cross-group comparison. | No |
| **NA** | No cross-group members (within-group-only group). | N/A |

### How cross-group classification works

1. **Codon fingerprinting** (`classify_sharing_status.py`): Compares
   consensus fingerprint of G1 members vs G2 members. Exact match +
   same family -> `ancestral`. Incompatible positions -> `independent`.
   Compatible position with conserved amino acid context -> `ancestral`.
   Compatible position but different protein context (legacy codon match
   only) -> `independent`. Different families -> `independent`.
   No fingerprints -> `uncertain`.

2. **Sequence refinement** (`compare_introner_sequences.py`): For groups
   still `uncertain`, uses cross-group body sequence identity as the only
   available evidence (>= 60% -> `likely_ancestral`, < 60% ->
   `likely_independent`). Also flags `ancestral` groups with identity
   < 60% as `ancestral_low_identity`.

### Why sequence identity doesn't override codon evidence

For cross-group comparisons with codon data, the codon-level classification
is treated as primary and sequence identity only flags outliers. Introners of
the same family share 60-80% body sequence identity whether or not they
descend from the same ancestral insertion. Two independent family-2
insertions in different genes will have similar sequences simply because they
are family-2 introners. Codon position and the flanking amino acid context
are the more direct evidence for the insertion event itself.

---

## Fingerprint matching details

Each introner's insertion site is characterized by a **hybrid key**:

1. **Flanking amino acid context**: 5 amino acids upstream + 5 downstream,
   translated from the CDS. Compared via sliding 6-mer window with 1
   mismatch tolerance.

2. **Legacy codon/intron position**: codon number + offset (CDS insertions)
   or intron number (intronic insertions). Fallback when aa context is
   unavailable or frameshifts cause different translations.

Two keys match if **either** the aa contexts match **or** the legacy
positions are compatible. This hybrid approach handles miniprot annotation
variability (different exon boundaries, non-multiple-of-3 CDS lengths).

For cross-group classification, an exact or sliding-window aa-context match
is required for an `ancestral` call. A legacy-only match (codon numbers
within tolerance but completely different protein context) is classified as
`independent` because the positional match is likely annotation noise.

---

## Current counts (pipeline run 2026-04-10)

### `within_group_status` (5953 ortholog groups)

| Category | Count | % |
|---|---|---|
| consistent | 5287 | 88.8% |
| singleton | 493 | 8.3% |
| low_identity | 164 | 2.8% |
| discordant | 7 | 0.1% |
| uncertain | 2 | <0.1% |

### `cross_group_status` (225 groups with G1+G2 members)

| Category | Count | % of cross-group |
|---|---|---|
| independent | 121 | 53.8% |
| ancestral | 56 | 24.9% |
| likely_ancestral | 27 | 12.0% |
| likely_independent | 15 | 6.7% |
| ancestral_low_identity | 6 | 2.7% |

5728 groups have `cross_group_status` = NA (within-group-only).

### Summary for downstream analysis

- **pi**: filter `within_group_status == 'consistent'` -> 5287 groups
- **dxy (ancestral/shared)**: `cross_group_status in ('ancestral', 'likely_ancestral')` -> 83 groups
  (+ 6 `ancestral_low_identity` for review)
- **dxy (independent)**: `cross_group_status in ('independent', 'likely_independent')` -> 136 groups
