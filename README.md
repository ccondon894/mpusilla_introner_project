# Micromonas pusilla Introner Project

Unified analysis pipeline for introner insertion pattern analysis in *Micromonas pusilla*.

## Overview

This project consolidates multiple analysis workflows for studying introner elements in *M. pusilla*:

- **Genotyping**: Introner presence/absence matrix generation
- **Evolution**: Diversity metrics and haplotype analysis
- **SNP Population Genetics**: Phylogenetics, demography, selection
- **Expression**: Isoform analysis, GLM modeling
- **GO Enrichment**: Functional annotation and enrichment analysis

## Quick Start

```bash
# Run full pipeline
cd /scratch1/chris/mpusilla_introner_project
snakemake --configfile config/config.yaml --cores 20

# Dry run to see what will be executed
snakemake -n --configfile config/config.yaml
```

## Directory Structure

```
├── config/           # Configuration files
├── data/             # Input data (symlinks)
├── workflow/         # Snakemake workflows
│   ├── Snakefile     # Master entry point
│   ├── rules/        # Individual rule files
│   └── envs/         # Conda environments
├── scripts/          # Python/R scripts
├── resources/        # Static resources
└── results/          # Output files
```

## Sample Groups

- **Group 1** (11 samples): CCMP1545 (reference), RCC114, RCC1614, RCC1698, RCC2482, RCC373, RCC465, RCC629, RCC692, RCC693, RCC833
- **Group 2** (2 samples): RCC1749 (reference), RCC3052

