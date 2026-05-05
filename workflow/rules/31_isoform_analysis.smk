# ============================================================
# 31_isoform_analysis.smk - Isoform Diversity Analysis
# ============================================================
#
# Analyzes isoform diversity and NMD (nonsense-mediated decay):
# 1. Shannon diversity analysis of isoform expression
# 2. NMD prediction analysis correlated with introner status
# 3. Isoform abundance characterization
#
# Adapted from:
# - /scratch1/chris/introner-expression-analysis/scripts/shannon_diversity_analysis_consolidated.py
# - /scratch1/chris/introner-expression-analysis/scripts/analyze_nmd_predictions.py
#
# ============================================================

import os
from pathlib import Path

# ============================================================
# CONFIGURATION
# ============================================================

# Output directories
EXPRESSION_DIR = RESULTS / "expression"
ISOFORM_DIR = EXPRESSION_DIR / "isoform_analysis"
DIVERSITY_DIR = ISOFORM_DIR / "diversity"
NMD_DIR = ISOFORM_DIR / "nmd"
COUNTS_DIR = EXPRESSION_DIR / "counts"
EXPRESSION_LOG_DIR = EXPRESSION_DIR / "logs"

# SQANTI3 output directory
SQANTI_DIR = EXPRESSION_DIR / "sqanti3"

# Mating type region to exclude
MT_SCAFFOLD = config["mating_type_region"]["scaffold"]
MT_START = config["mating_type_region"]["start"]
MT_END = config["mating_type_region"]["end"]


# ============================================================
# SHANNON DIVERSITY ANALYSIS
# ============================================================

rule calculate_shannon_diversity:
    """
    Calculate Shannon diversity index for isoform expression.

    For each gene, calculates:
    - Shannon Index (H): H = -sum(p_i * ln(p_i))
    - Shannon Evenness (J'): J' = H / ln(k)

    Higher values indicate more even isoform expression.
    Filters out mating-type region genes.
    """
    input:
        sqanti_data = SQANTI_DIR / "parsed_sqanti3_data.tsv",
        gtf = ANNOTATIONS_DIR / "CCMP1545.gtf"
    output:
        diversity = DIVERSITY_DIR / "shannon_diversity_per_gene.csv",
        summary = DIVERSITY_DIR / "shannon_diversity_summary.txt",
        plot = FIGURES_DIR / "shannon_diversity_distribution.pdf"
    params:
        mt_scaffold = MT_SCAFFOLD,
        mt_start = MT_START,
        mt_end = MT_END
    log:
        EXPRESSION_LOG_DIR / "shannon_diversity.log"
    shell:
        """
        mkdir -p {DIVERSITY_DIR}

        python {PROJECT_ROOT}/scripts/expression/shannon_diversity_analysis_consolidated.py \
            --sqanti {input.sqanti_data} \
            --gtf {input.gtf} \
            --output {output.diversity} \
            --summary {output.summary} \
            --plot {output.plot} \
            --mt_scaffold {params.mt_scaffold} \
            --mt_start {params.mt_start} \
            --mt_end {params.mt_end} \
            2> {log}
        """


rule diversity_by_introner_status:
    """
    Compare Shannon diversity between genes with/without introners.

    Tests whether introner presence correlates with isoform diversity.
    """
    input:
        diversity = DIVERSITY_DIR / "shannon_diversity_per_gene.csv",
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv"
    output:
        comparison = DIVERSITY_DIR / "diversity_by_introner_status.csv",
        plot = FIGURES_DIR / "diversity_introner_comparison.pdf",
        stats = DIVERSITY_DIR / "diversity_introner_stats.txt"
    log:
        EXPRESSION_LOG_DIR / "diversity_by_introner.log"
    shell:
        """
        python {PROJECT_ROOT}/scripts/expression/shannon_diversity_analysis_consolidated.py \
            --mode compare_introner \
            --diversity {input.diversity} \
            --genotype_matrix {input.genotype_matrix} \
            --output {output.comparison} \
            --plot {output.plot} \
            --stats {output.stats} \
            2> {log}
        """


# ============================================================
# INTRONER LOCI BED FILE
# ============================================================

rule create_introner_loci_bed:
    """
    Extract CCMP1545 introner loci from genotype matrix into BED format.
    """
    input:
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv"
    output:
        bed = GENOTYPING_DIR / "introner_loci.bed"
    run:
        import pandas as pd
        df = pd.read_csv(input.matrix, sep='\t')
        ref = df[(df['sample'] == 'CCMP1545') & (df['presence'] == 1)]
        bed = ref[['contig', 'start', 'end', 'gene', 'ortholog_id', 'family']].copy()
        bed['start'] = bed['start'].astype(int)
        bed['end'] = bed['end'].astype(int)
        bed['family'] = bed['family'].astype(int)
        bed.to_csv(output.bed, sep='\t', header=False, index=False)


# ============================================================
# NMD (NONSENSE-MEDIATED DECAY) ANALYSIS
# ============================================================

rule analyze_nmd_predictions:
    """
    Analyze NMD predictions from SQANTI3 classification.

    Correlates NMD status with introner presence to test whether
    introner loss allows more NMD-targeted isoforms to survive.

    Uses Fisher's exact test for statistical significance.
    """
    input:
        sqanti_data = SQANTI_DIR / "parsed_sqanti3_data.tsv",
        gtf = ANNOTATIONS_DIR / "CCMP1545.gtf",
        introner_bed = GENOTYPING_DIR / "introner_loci.bed"
    output:
        analysis = NMD_DIR / "nmd_introner_analysis.txt",
        contingency = NMD_DIR / "nmd_contingency_tables.csv",
        plot = FIGURES_DIR / "nmd_introner_association.pdf"
    params:
        # Include only high-confidence structural categories
        categories = "full-splice_match,novel_in_catalog,novel_not_in_catalog"
    log:
        EXPRESSION_LOG_DIR / "nmd_analysis.log"
    shell:
        """
        mkdir -p {NMD_DIR}

        python {PROJECT_ROOT}/scripts/expression/analyze_nmd_predictions.py \
            --sqanti {input.sqanti_data} \
            --gtf {input.gtf} \
            --introner_bed {input.introner_bed} \
            --output {output.analysis} \
            --contingency {output.contingency} \
            --plot {output.plot} \
            --categories {params.categories} \
            2> {log}
        """


rule nmd_by_strain:
    """
    Analyze NMD predictions stratified by strain.

    Compares NMD patterns across CCMP1545, RCC1614, and RCC1749.
    """
    input:
        sqanti_data = SQANTI_DIR / "parsed_sqanti3_data.tsv",
        gtf_ccmp1545 = ANNOTATIONS_DIR / "CCMP1545.gtf",
        gtf_rcc1614 = ANNOTATIONS_DIR / "RCC1614.gtf",
        gtf_rcc1749 = get_gtf("RCC1749")
    output:
        analysis = NMD_DIR / "nmd_by_strain_analysis.txt",
        plot = FIGURES_DIR / "nmd_by_strain.pdf"
    log:
        EXPRESSION_LOG_DIR / "nmd_by_strain.log"
    shell:
        """
        python {PROJECT_ROOT}/scripts/expression/analyze_nmd_predictions.py \
            --mode by_strain \
            --sqanti {input.sqanti_data} \
            --gtf_ccmp1545 {input.gtf_ccmp1545} \
            --gtf_rcc1614 {input.gtf_rcc1614} \
            --gtf_rcc1749 {input.gtf_rcc1749} \
            --output {output.analysis} \
            --plot {output.plot} \
            2> {log}
        """


# ============================================================
# ISOFORM COUNT ANALYSIS
# ============================================================

rule count_isoforms_per_gene:
    """
    Count number of isoforms per gene for each strain.
    """
    input:
        sqanti_data = SQANTI_DIR / "parsed_sqanti3_data.tsv"
    output:
        counts = ISOFORM_DIR / "isoform_counts_per_gene.csv"
    log:
        EXPRESSION_LOG_DIR / "count_isoforms.log"
    shell:
        """
        mkdir -p {ISOFORM_DIR}

        python {PROJECT_ROOT}/scripts/expression/isoform_analysis/isoform_abundance_analysis.py \
            --sqanti {input.sqanti_data} \
            --output {output.counts} \
            --mode count \
            2> {log} || \
        python -c "
import pandas as pd
# Fallback: simple isoform counting
df = pd.read_csv('{input.sqanti_data}', sep='\\t')
counts = df.groupby(['gene_id', 'strain']).size().reset_index(name='n_isoforms')
counts.to_csv('{output.counts}', index=False)
"
        """


# ============================================================
# DATA PREPARATION FOR GLM MODELING
# ============================================================

rule prepare_isoform_data:
    """
    Prepare comprehensive isoform data for GLM modeling.

    Merges:
    - Isoform counts per gene
    - Introner gain/loss features (from genotype matrix)
    - CDS lengths (from GTF)
    - Mean expression (from featureCounts counts matrix)
    """
    input:
        isoform_counts = ISOFORM_DIR / "isoform_counts_per_gene.csv",
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        gtf = ANNOTATIONS_DIR / "CCMP1545.gtf",
        counts_matrix = COUNTS_DIR / "merged_counts_matrix.csv"
    output:
        data = ISOFORM_DIR / "isoform_introner_data.csv",
        filtered = ISOFORM_DIR / "isoform_introner_data_filtered.csv"
    params:
        mt_scaffold = MT_SCAFFOLD,
        mt_start = MT_START,
        mt_end = MT_END
    log:
        EXPRESSION_LOG_DIR / "prepare_isoform_data.log"
    shell:
        """
        python {PROJECT_ROOT}/scripts/expression/prepare_isoform_data.py \
            --isoform_counts {input.isoform_counts} \
            --genotype_matrix {input.genotype_matrix} \
            --gtf {input.gtf} \
            --output {output.data} \
            --filtered {output.filtered} \
            --counts_matrix {input.counts_matrix} \
            --mt_scaffold {params.mt_scaffold} \
            --mt_start {params.mt_start} \
            --mt_end {params.mt_end} \
            2> {log}
        """


# ============================================================
# TARGET RULES
# ============================================================

rule isoform_analysis_complete:
    """
    Target: Complete isoform diversity and NMD analysis.
    """
    input:
        DIVERSITY_DIR / "shannon_diversity_per_gene.csv",
        DIVERSITY_DIR / "diversity_by_introner_status.csv",
        NMD_DIR / "nmd_introner_analysis.txt",
        ISOFORM_DIR / "isoform_introner_data_filtered.csv",
        FIGURES_DIR / "shannon_diversity_distribution.pdf",
        FIGURES_DIR / "nmd_introner_association.pdf"


rule shannon_diversity_only:
    """
    Target: Shannon diversity analysis only.
    """
    input:
        DIVERSITY_DIR / "shannon_diversity_per_gene.csv",
        FIGURES_DIR / "shannon_diversity_distribution.pdf"


rule nmd_analysis_only:
    """
    Target: NMD analysis only.
    """
    input:
        NMD_DIR / "nmd_introner_analysis.txt",
        FIGURES_DIR / "nmd_introner_association.pdf"
