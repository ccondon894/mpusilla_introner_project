# PRD: Micromonas pusilla Introner Research Project Refactoring

**Project:** mpusilla_introner_project
**Author:** Chris
**Created:** 2026-01-23
**Status:** Planning
**Location:** `/scratch1/chris/mpusilla_introner_project/`

---

## Executive Summary

Refactor the *Micromonas pusilla* introner research project from 8+ scattered directories into a unified, reproducible project with full Snakemake automation.

### Goals
- Consolidate all analyses into a single project directory
- Eliminate circular dependencies (specifically the GTF introner annotation step)
- Convert manual shell scripts to Snakemake workflows (especially SNP popgen)
- Preserve all original code (read-only references)

### Non-Goals
- Re-running computationally expensive steps (use existing outputs where valid)
- Including deprecated workflows (ortholog_detection_pipeline v1)
- Including dropped samples (RCC835, RCC647, CCMP490)

---

## Progress Tracking

A `progress.txt` file will be maintained at:
`/scratch1/chris/mpusilla_introner_project/progress.txt`

Format:
```
# Progress Tracker - mpusilla_introner_project
# Updated: YYYY-MM-DD

## Iteration 1: Project Foundation & Core Genotyping
- [ ] 1.1 Directory structure created
- [ ] 1.2 config.yaml created
- [ ] 1.3 sample_metadata.tsv created
- [ ] 1.4 Data symlinks established
- [ ] 1.5 resources copied
- [ ] 2.1 00_data_prep.smk created
- [ ] 2.2 03_gtf_liftover.smk created
- [ ] 2.3 04_blast_genotyping.smk created
- [ ] 2.4 05_ortholog_detection.smk created
- [ ] 2.5 06_coverage_calling.smk created
- [ ] 2.6 Scripts copied to scripts/genotyping/
- [ ] 2.7 Master Snakefile created
- [ ] 2.8 VERIFIED: Dry run successful

## Iteration 2: Evolution Analysis
- [ ] 3.1 12_genome_alignment_circos.smk created
- [ ] 3.2 Scripts copied to scripts/genome_alignment/
- [ ] 3.3 10_evolution_analysis.smk created
- [ ] 3.4 11_haplotype_diversity.smk created
- [ ] 3.5 Scripts copied to scripts/evolution/
- [ ] 3.6 VERIFIED: Dry run successful

## Iteration 3: SNP Population Genetics
- [ ] 4.1 20_snp_calling.smk created
- [ ] 4.2 21_snp_popgen.smk created
- [ ] 4.3 22_phylogenetics.smk created
- [ ] 4.4 23_demography.smk created
- [ ] 4.5 24_selection.smk created
- [ ] 4.6 Scripts copied to scripts/popgen/
- [ ] 4.7 VERIFIED: Dry run successful

## Iteration 4: Expression & GO Analysis
- [ ] 5.1 30_expression.smk created
- [ ] 5.2 31_isoform_analysis.smk created
- [ ] 5.3 32_glm_modeling.smk created
- [ ] 5.4 Scripts copied to scripts/expression/
- [ ] 6.1 40_go_enrichment.smk created
- [ ] 6.2 Scripts copied to scripts/go_analysis/
- [ ] 6.3 VERIFIED: Full pipeline dry run successful
```

---

## Directory Structure

```
/scratch1/chris/mpusilla_introner_project/
├── PRD.md                             # This document
├── progress.txt                       # Implementation progress tracker
├── README.md                          # Project overview for users
│
├── config/
│   └── config.yaml                    # Central configuration (paths, params, samples)
│
├── data/                              # Input data (symlinks to originals)
│   ├── reads/                         # → /scratch1/alex/pusilla/data/reads/
│   ├── rna_seq/                       # → /storage1/chris/introner_RNA_seq/fastq
│   ├── references/
│   │   ├── CCMP1545_v3.fa             # → /scratch1/alex/pusilla/data/ref/v3/...
│   │   └── RCC1749_assembly.fa        # → /scratch1/alex/pusilla/data/intronerless/...
│   ├── external_db/                   # GO, PANTHER, InterPro files
│   └── r2c2_bams/                     # → /scratch2/russ/introner/...
│
├── workflow/
│   ├── Snakefile                      # Master workflow entry point
│   ├── rules/
│   │   ├── 00_data_prep.smk
│   │   ├── 03_gtf_liftover.smk
│   │   ├── 04_blast_genotyping.smk
│   │   ├── 05_ortholog_detection.smk
│   │   ├── 06_coverage_calling.smk
│   │   ├── 10_evolution_analysis.smk
│   │   ├── 11_haplotype_diversity.smk
│   │   ├── 12_genome_alignment_circos.smk
│   │   ├── 20_snp_calling.smk
│   │   ├── 21_snp_popgen.smk
│   │   ├── 22_phylogenetics.smk
│   │   ├── 23_demography.smk
│   │   ├── 24_selection.smk
│   │   ├── 30_expression.smk
│   │   ├── 31_isoform_analysis.smk
│   │   ├── 32_glm_modeling.smk
│   │   └── 40_go_enrichment.smk
│   └── envs/
│       ├── genotyping.yaml
│       ├── popgen.yaml
│       └── expression.yaml
│
├── scripts/
│   ├── genotyping/
│   ├── evolution/
│   ├── genome_alignment/
│   ├── popgen/
│   ├── expression/
│   └── go_analysis/
│
├── resources/
│   ├── introner_truth_set.fa
│   ├── go.obo
│   └── sample_metadata.tsv
│
└── results/
    ├── assemblies/
    ├── annotations/
    ├── genotyping/
    ├── evolution/
    ├── genome_alignment/
    ├── snp_popgen/
    ├── expression/
    ├── go_enrichment/
    └── figures/
```

---

## Iteration 1: Project Foundation & Core Genotyping

### 1.1 Project Skeleton

**Objective:** Create directory structure and configuration files.

**Tasks:**

| Task | Description | Output |
|------|-------------|--------|
| 1.1.1 | Create directory tree | All directories exist |
| 1.1.2 | Create config.yaml | `config/config.yaml` |
| 1.1.3 | Create sample_metadata.tsv | `resources/sample_metadata.tsv` |
| 1.1.4 | Create symlinks to input data | `data/` populated |
| 1.1.5 | Copy introner_truth_set.fa | `resources/introner_truth_set.fa` |
| 1.1.6 | Copy go.obo | `resources/go.obo` |

**config.yaml structure:**
```yaml
# Sample configuration
samples:
  group1: [CCMP1545, RCC114, RCC1614, RCC1698, RCC2482, RCC373, RCC465, RCC629, RCC692, RCC693, RCC833]
  group2: [RCC1749, RCC3052]
  reference: CCMP1545

# Input paths (symlinked)
paths:
  reads: data/reads
  rna_seq: data/rna_seq
  references:
    ccmp1545: data/references/CCMP1545_v3.fa
    rcc1749: data/references/RCC1749_assembly.fa

# Tool parameters
params:
  blast:
    identity_threshold: 0.8
    similarity_cutoff: 0.7
  bwa:
    mapq_threshold: 30
    match_length: 70
  bwa_rescue:
    mapq_threshold: 20
    match_length: 50
```

**sample_metadata.tsv:**
```
sample_id	group	is_reference	notes
CCMP1545	group1	true	Reference strain (834)
RCC114	group1	false
RCC1614	group1	false	Has RNA-seq
RCC1749	group2	true	Group2 reference, has RNA-seq
RCC3052	group2	false
...
```

### 1.2 Core Genotyping Pipeline

**Objective:** Consolidate genotyping workflows (rules 00-06).

**Source Files:**
- `/scratch1/chris/introner-genotyping-pipeline/rules/blast_genotyping_graph.smk`
- `/scratch1/chris/introner-genotyping-pipeline/rules/ortholog_detection_pipeline_v2.smk`
- `/scratch1/chris/introner-genotyping-pipeline/rules/call_introner_presence.smk`
- `/scratch1/chris/introner-expression-analysis/rules/make_sample_gtfs.smk`

**Rules to Create:**

#### 00_data_prep.smk
- Index reference genomes (BWA, samtools faidx)
- Create genome dictionaries
- Validate input files exist

#### 03_gtf_liftover.smk
Adapted from `make_sample_gtfs.smk`:
- Extract reference proteins
- Run miniprot alignment
- Resolve paralog mappings
- Filter GTF
- Rename gene IDs
- **EXCLUDED:** `update_gtf_with_introners.py` (removes circular dependency)

#### 04_blast_genotyping.smk
Adapted from `blast_genotyping_graph.smk`:
- BLAST introner truth set against sample genomes
- Filter by identity, length, complexity
- Extract candidate loci with flanks
- Validate against reference annotations

#### 05_ortholog_detection.smk
Adapted from `ortholog_detection_pipeline_v2.smk`:
- Phase 0: GTF to BED conversion
- Phase 1A: Synteny mapping (upstream/downstream genes)
- Phase 1B: Flank extraction, BWA indexing
- Phase 1C: Initial ortholog pairing
- Phase 2: Rescue alignments for Group1↔Group2 pairs
- Build genotype matrix
- Fix orientations
- Annotate missing data

#### 06_coverage_calling.smk
Adapted from `call_introner_presence.smk`:
- Align reads to sample genomes
- Calculate coverage at introner loci
- Update genotype matrix with coverage calls

**Scripts to Copy:**
```
Source: /scratch1/chris/introner-genotyping-pipeline/scripts/
Target: scripts/genotyping/

- introner_genotype_matrix_builder.py
- fix_matrix_orientations.py
- annotate_missing_data.py
- validate_candidate_loci.py
- extract_candidate_flanks.py
```

```
Source: /scratch1/chris/introner-expression-analysis/scripts/
Target: scripts/genotyping/

- resolve_paralog_mappings.py
- filter_protein_gtf.py
- rename_gene_ids.py
- gtf_stats.py
```

**Verification:**
```bash
cd /scratch1/chris/mpusilla_introner_project
snakemake -s workflow/Snakefile --configfile config/config.yaml -n --until genotype_matrix
```

---

## Iteration 2: Evolution Analysis

### 2.1 Whole Genome Alignment & Circos Visualization

**Objective:** Port the MUMmer whole-genome alignment and Circos synteny visualization comparing CCMP1545 vs RCC1749.

**Source File:** `/scratch1/chris/introner-group-comparison-analysis/rules/mummer_circos.smk`

#### 12_genome_alignment_circos.smk

Rules (in dependency order):

1. **nucmer** — `nucmer` alignment of CCMP1545 vs RCC1749 vg_paths assemblies, filtered with `delta-filter -1`
2. **mummerplot** — dot plot PNG from `.1delta` file
3. **generate_karyotype** — build `.kar` from both `.fa.fai` indices
4. **generate_links** — extract synteny blocks (min-length 1000 bp, min-identity 80%)
5. **analyze_synteny** — reorder karyotype by synteny
6. **flip_chromosomes** — flip RCC1749 contig orientation
7. **filter_reorder** — remove contigs `2,1,33,30,32,40,43,24,42`; move contig 20 after 12 before 25
8. **generate_ribbons** — colored synteny ribbons (gap_threshold=50000, merge_threshold=750)
9. **generate_circos_config** — produce `circos.conf` dynamically via `generate_circos_config.py`
10. **render_circos** — `circos -conf circos.conf` → `circos.png`, `circos.svg`

**Path Mapping (old → new):**

| Old | New |
|-----|-----|
| `ASSEMBLY_DIR/{g}.vg_paths.fa` | `ASSEMBLIES_DIR / "{g}.vg_paths.fa"` (from `copy_assemblies`) |
| `ASSEMBLY_DIR/{g}.vg_paths.fa.fai` | `ASSEMBLIES_DIR / "{g}.vg_paths.fa.fai"` (from `index_assembly_fasta`) |
| `MUMMER_DIR` | `GENOME_ALIGNMENT_DIR / "mummer"` |
| `FIG_DIR` (dot plot) | `FIGURES_DIR / "genome_alignment"` |
| `CIRCOS_DIR` | `GENOME_ALIGNMENT_DIR / "circos"` |
| `scripts/...` | `scripts/genome_alignment/...` |

**New Snakefile additions:**
- Add `GENOME_ALIGNMENT_DIR = RESULTS / "genome_alignment"` to directory variables
- Include `rules/12_genome_alignment_circos.smk` inside the `LOAD_EVOLUTION` block
- Add `genome_alignment_circos` target rule inside the `if LOAD_EVOLUTION:` block

**config.yaml additions** (under `params:`):
```yaml
mummer:
  min_link_length: 1000
  min_identity: 80.0
circos:
  gap_threshold: 50000
  merge_threshold: 750
  remove_contigs: "2,1,33,30,32,40,43,24,42"
  move_contig: 20
  insert_after: 12
  insert_before: 25
```

**Scripts to Copy:**
```
Source: /scratch1/chris/introner-group-comparison-analysis/scripts/
Target: scripts/genome_alignment/

- generate_karyotype.py
- generate_links.py
- analyze_synteny.py
- flip_rcc1749_order.py
- filter_reorder_karyotype.py
- generate_ribbon_links_with_colors.py
- generate_circos_config.py
```

**Key outputs:**
- `results/genome_alignment/mummer/CCMP1545_vs_RCC1749.1delta`
- `results/figures/genome_alignment/CCMP1545_vs_RCC1749.png`
- `results/genome_alignment/circos/circos.png`
- `results/genome_alignment/circos/circos.svg`

**Verification:**
```bash
snakemake --configfile config/config.yaml --config modules=evolution -n genome_alignment_circos
snakemake --configfile config/config.yaml --config modules=evolution --cores 8 genome_alignment_circos
```

---

### 2.2 Evolution & Haplotype Workflows

**Objective:** Consolidate evolution and haplotype diversity analyses.

**Source Files:**
- `/scratch1/chris/introner-genotyping-pipeline/rules/all_samples_evolution_analysis.smk`
- `/scratch1/chris/introner-genotyping-pipeline/rules/group1_evolution_analysis.smk`
- `/scratch1/chris/introner-genotyping-pipeline/rules/haplotype_diversity_workflow.smk`

**Rules to Create:**

#### 10_evolution_analysis.smk
- Build consensus sequences for flanking regions
- Classify orthologs by fixation category
- Run MAFFT alignments
- Calculate diversity metrics (π, dS)
- Generate visualizations

#### 11_haplotype_diversity.smk
- Parse 4-fold degenerate site VCF
- Calculate PHDR/AHDR ratios
- Generate frequency-based plots

**Scripts to Copy:**
```
Source: /scratch1/chris/introner-genotyping-pipeline/scripts/
Target: scripts/evolution/

- build_consensus_sequences.py
- classify_fixation_categories.py
- calculate_diversity_metrics.py
- plot_diversity_boxplots.py
```

**Verification:**
```bash
snakemake -n --until haplotype_diversity_ratios
```

---

## Iteration 3: SNP Population Genetics (NEW AUTOMATION)

### 3.1 Overview

**Objective:** Convert manual `commands.sh` scripts to Snakemake workflows.

**Source Directory:** `/scratch1/chris/introner_vis/`

This is the largest new automation effort. Currently these analyses are run via shell scripts with manual intervention.

### 3.2 Rules to Create

#### 20_snp_calling.smk
**Source:** `/scratch1/chris/introner-genotyping-final/rules/call_haplotypes_to_CCMP1545.smk`

Using simplified direct alignment approach:
1. Index CCMP1545 reference
2. Align FASTQs with BWA
3. Add read groups (GATK)
4. Run HaplotypeCaller (haploid mode)
5. Combine GVCFs
6. Joint genotyping

**Output:** `results/snp_popgen/vcf/mpusilla.joint.vcf.gz`

#### 21_snp_popgen.smk
**Source:** `/scratch1/chris/introner_vis/data_prep/commands.sh`, `/scratch1/chris/introner_vis/degenotate/`

1. Run degenotate for 4-fold sites
2. Extract SNPs from VCF
3. Filter by 4-fold degeneracy
4. Create sample subsets (intronerful, no missing data)
5. Run SnpEff annotation

**Outputs:**
- `results/snp_popgen/vcf/mpusilla.snps.4d.notMT.vcf.gz`
- `results/snp_popgen/vcf/mpusilla.snps.4d.notMT.intronerful.vcf.gz`
- `results/snp_popgen/vcf/mpusilla.snps.snpEff.vcf`

#### 22_phylogenetics.smk
**Source:** `/scratch1/chris/introner_vis/iqtree/`

1. Convert VCF to PHYLIP (`vcf2phylip.py`)
2. Run IQ-TREE ML
3. Reroot tree
4. Plot tree (`plot_tree_improved.py`)

**Output:** `results/snp_popgen/phylogenetics/phylogenetic_tree.pdf` (PUBLICATION FIGURE)

#### 23_demography.smk
**Source:** `/scratch1/chris/introner_vis/dadi/`

1. Convert VCF to dadi format
2. Fit demographic models (`fit_model.py`)
3. Visualize demography (`visualize_demography.py`)

**Output:** `results/snp_popgen/demography/demographic_model.pdf` (MAIN TEXT FIGURE)

#### 24_selection.smk
**Source:** `/scratch1/chris/introner_vis/basic_popgen/`, `/scratch1/chris/introner_vis/pca/`, `/scratch1/chris/introner_vis/ld/`, `/scratch1/chris/introner_vis/polymorphism_analysis/`, `/scratch1/chris/introner_vis/recombination_analysis/`

1. **Basic popgen:** Tajima's D, 2D AFS (MAIN FIGURE), 1D AFS
2. **PCA:** plink LD pruning, eigendecomposition, plotting
3. **LD:** Compute and plot linkage disequilibrium
4. **SFS:** Unfolded SFS density plots (MAIN FIGURE)
5. **Recombination:** pyrho comparison plots (MAIN FIGURE)

**Scripts to Copy:**
```
Source: /scratch1/chris/introner_vis/
Target: scripts/popgen/

basic_popgen/
- tajima_d.py
- 2D_afs.py
- unfolded_1D_afs.py

dadi/
- fit_model.py
- visualize_demography.py
- plot_dadi.py

iqtree/
- vcf2phylip.py
- plot_tree_improved.py

ld/
- compute_ld.py
- summarize_ld.py
- plot_ld.py

pca/
- plot_PCA.py
- plot_PCA.r

polymorphism_analysis/
- phase2_sfs_analysis.py

recombination_analysis/
- compare_recombination_gene_exonic_introners.py
- compare_recombination_gene_exonic_loss_gain.py
- compare_recombination_gene_exonic_frequency_based.py
- compare_recombination_gene_exonic_group2_introners.py
- plot_recombination_boxplots.py
```

**Verification:**
```bash
snakemake -n --until snp_popgen_complete
```

---

## Iteration 4: Expression & GO Analysis

### 4.1 Expression Analysis

**Objective:** Consolidate expression and isoform analyses.

**Source Directory:** `/scratch1/chris/introner-expression-analysis/`

#### 30_expression.smk
- Run SQANTI3 on R2C2 data
- Run featureCounts for quantification

#### 31_isoform_analysis.smk
- Shannon diversity analysis (`shannon_diversity_analysis_consolidated.py`)
- NMD analysis (`analyze_nmd_predictions.py`)

#### 32_glm_modeling.smk
- Poisson GLM for isoform diversity
- Negative binomial GLM for expression levels

**Scripts to Copy:**
```
Source: /scratch1/chris/introner-expression-analysis/
Target: scripts/expression/

scripts/
- shannon_diversity_analysis_consolidated.py
- analyze_nmd_predictions.py

isoform_analysis/poisson_modeling/
- poisson_GLM_regression.py
- negative_binomial_glm.py
```

### 4.2 GO Enrichment

**Objective:** Automate GO annotation enhancement and enrichment analysis.

**Source Directory:** `/scratch1/chris/mpusilla_go_analysis/`

#### 40_go_enrichment.smk

Steps:
1. Generate GO mapping files
   - pfam_to_go.py
   - ko_to_go.py (requires internet)
   - tair_to_go.py
   - panther_to_go.py
2. Enhance GO coverage (`enhance_go_coverage.py`)
3. Phase 1 gene classification (`introner_phase1_analysis.py`)
4. Phase 2 enrichment (`introner_phase2_enrichment.py`)
5. Phase 2 family enrichment (`introner_phase2_family_enrichment.py`)

**Scripts to Copy:**
```
Source: /scratch1/chris/mpusilla_go_analysis/scripts/
Target: scripts/go_analysis/

- pfam_to_go.py
- ko_to_go.py
- tair_to_go.py
- panther_to_go.py
- enhance_go_coverage.py
- introner_phase1_analysis.py
- introner_phase2_enrichment.py
- introner_phase2_family_enrichment.py
```

**External Data Needed:**
- `data/external_db/interpro2go`
- `data/external_db/gene_association.tair`
- `data/external_db/PANTHER19.0_HMM_classifications`

**Verification:**
```bash
snakemake -n --until go_enrichment_complete
```

---

## Dependency Graph

```
                         ┌─────────────────────┐
                         │    INPUT DATA       │
                         │ (reads, refs, R2C2) │
                         └──────────┬──────────┘
                                    │
           ┌────────────────────────┼────────────────────────┐
           │                        │                        │
           ▼                        ▼                        ▼
    ┌─────────────┐         ┌─────────────┐         ┌─────────────┐
    │ 03_gtf_     │         │ 04_blast_   │         │ 20_snp_     │
    │ liftover    │         │ genotyping  │         │ calling     │
    └──────┬──────┘         └──────┬──────┘         └──────┬──────┘
           │                        │                        │
           └────────────┬───────────┘                        │
                        ▼                                    ▼
               ┌─────────────────┐                  ┌─────────────────┐
               │ 05_ortholog_    │                  │ 21_snp_popgen   │
               │ detection       │                  │ 22_phylogenetics│
               └────────┬────────┘                  │ 23_demography   │
                        │                           │ 24_selection    │
                        ▼                           └─────────────────┘
               ┌─────────────────┐
               │ 06_coverage_    │
               │ calling         │
               └────────┬────────┘
                        │
                        ▼
              ┌─────────────────────┐
              │   GENOTYPE MATRIX   │
              └──────────┬──────────┘
                         │
        ┌────────────────┼────────────────┐
        │                │                │
        ▼                ▼                ▼
 ┌─────────────┐  ┌─────────────┐  ┌─────────────┐
 │ 10_evolution│  │ 30_express. │  │ 40_go_      │
 │ 11_haplotype│  │ 31_isoform  │  │ enrichment  │
 └─────────────┘  │ 32_glm      │  └─────────────┘
                  └─────────────┘
```

---

## Excluded Items

### Deprecated Workflows
- `ortholog_detection_pipeline.smk` (v1) - use v2 only
- `group1_4fold_analysis.smk` - superseded

### Excluded Samples
- RCC835, RCC647, CCMP490 (dropped from analysis)

### Removed Steps
- `update_gtf_with_introners.py` (was causing circular dependency)

### Not Migrated
- Old dated genotype matrix files
- Temporary/exploratory directories
- Unused scripts

---

## Verification Checklist

After each iteration, verify:

1. **Dry run succeeds:**
   ```bash
   snakemake -n --configfile config/config.yaml
   ```

2. **DAG is acyclic:**
   ```bash
   snakemake --dag | dot -Tpng > dag.png
   ```

3. **Outputs match originals (for key files):**
   - `genotype_matrix.tsv` vs `genotype_matrix_updated_12112025.tsv`
   - Publication figures vs existing versions

4. **Update progress.txt** with completed items

---

## Risk Assessment

| Risk | Mitigation |
|------|------------|
| Missing dependencies | Document all external tools in conda envs |
| Path hardcoding | Use config.yaml for all paths |
| Script incompatibility | Test each script independently before integration |
| Large compute time | Use existing outputs where appropriate |

---

## Appendix: Source Directory Mapping

| New Location | Original Source |
|--------------|-----------------|
| `workflow/rules/03_gtf_liftover.smk` | `/scratch1/chris/introner-expression-analysis/rules/make_sample_gtfs.smk` |
| `workflow/rules/04_blast_genotyping.smk` | `/scratch1/chris/introner-genotyping-pipeline/rules/blast_genotyping_graph.smk` |
| `workflow/rules/05_ortholog_detection.smk` | `/scratch1/chris/introner-genotyping-pipeline/rules/ortholog_detection_pipeline_v2.smk` |
| `workflow/rules/06_coverage_calling.smk` | `/scratch1/chris/introner-genotyping-pipeline/rules/call_introner_presence.smk` |
| `workflow/rules/10_evolution_analysis.smk` | `/scratch1/chris/introner-genotyping-pipeline/rules/all_samples_evolution_analysis.smk` + `group1_evolution_analysis.smk` |
| `workflow/rules/11_haplotype_diversity.smk` | `/scratch1/chris/introner-genotyping-pipeline/rules/haplotype_diversity_workflow.smk` |
| `workflow/rules/12_genome_alignment_circos.smk` | `/scratch1/chris/introner-group-comparison-analysis/rules/mummer_circos.smk` |
| `workflow/rules/20_snp_calling.smk` | `/scratch1/chris/introner-genotyping-final/rules/call_haplotypes_to_CCMP1545.smk` |
| `workflow/rules/21-24_*.smk` | `/scratch1/chris/introner_vis/` (NEW automation) |
| `workflow/rules/30-32_*.smk` | `/scratch1/chris/introner-expression-analysis/` |
| `workflow/rules/40_go_enrichment.smk` | `/scratch1/chris/mpusilla_go_analysis/` |
