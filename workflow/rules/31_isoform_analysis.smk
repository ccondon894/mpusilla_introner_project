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
ISOFORM_TURNOVER_DIR = ISOFORM_DIR / "turnover"
DIVERSITY_DIR = ISOFORM_DIR / "diversity"
NMD_DIR = ISOFORM_DIR / "nmd"
COUNTS_DIR = EXPRESSION_DIR / "counts"
EXPRESSION_LOG_DIR = EXPRESSION_DIR / "logs"

# SQANTI3 output directory
SQANTI_DIR = EXPRESSION_DIR / "sqanti3"

# Representative genes for the introner gene/transcript-structure track figure.
ISOFORM_TRACK_GENES = config.get(
    "isoform_track_genes",
    [
        "MicpuC2.est_orfs.13_5950_4275896:1.3.0.228",
        "estExt_Genewise1Plus.C_3_t40031.3.0.228",
    ],
)
ISOFORM_TRACK_PDFS = expand(
    FIGURES_DIR / "introner_isoform_tracks" / "{gene}.isoform_model_tracks.pdf",
    gene=ISOFORM_TRACK_GENES,
)
ISOFORM_TRACK_PNGS = expand(
    FIGURES_DIR / "introner_isoform_tracks" / "{gene}.isoform_model_tracks.png",
    gene=ISOFORM_TRACK_GENES,
)

# Mating type region to exclude
MT_SCAFFOLD = config["mating_type_region"]["scaffold"]
MT_START = config["mating_type_region"]["start"]
MT_END = config["mating_type_region"]["end"]
TURNOVER_MT_CONTIG = config["mating_type_region"]["group1"]["contig"]


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
        plot = FIGURES_DIR / "shannon_diversity_distribution.pdf",
        png = FIGURES_DIR / "shannon_diversity_distribution.png"
    params:
        mt_scaffold = MT_SCAFFOLD,
        mt_start = MT_START,
        mt_end = MT_END
    log:
        EXPRESSION_LOG_DIR / "shannon_diversity.log"
    conda: "../envs/isoform_analysis.yaml"
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
        png = FIGURES_DIR / "diversity_introner_comparison.png",
        stats = DIVERSITY_DIR / "diversity_introner_stats.txt"
    log:
        EXPRESSION_LOG_DIR / "diversity_by_introner.log"
    conda: "../envs/isoform_analysis.yaml"
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


rule create_strain_introner_loci_bed:
    """
    Extract strain-specific introner loci from the genotype matrix.

    These BEDs are used for per-strain NMD analysis so each SQANTI3
    transcript set is compared against introners in its own genome frame.
    """
    input:
        matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv"
    output:
        bed = GENOTYPING_DIR / "introner_loci" / "{sample}.introner_loci.bed"
    run:
        import pandas as pd
        df = pd.read_csv(input.matrix, sep='\t')
        ref = df[(df['sample'] == wildcards.sample) & (df['presence'] == 1)].copy()
        bed = ref[['contig', 'start', 'end', 'gene', 'ortholog_id', 'family']].copy()
        bed = bed.dropna(subset=['contig', 'start', 'end'])
        bed['start'] = bed['start'].astype(int)
        bed['end'] = bed['end'].astype(int)
        bed['family'] = bed['family'].fillna(-1).astype(int)
        bed.to_csv(output.bed, sep='\t', header=False, index=False)


# ============================================================
# NMD (NONSENSE-MEDIATED DECAY) ANALYSIS
# ============================================================

rule analyze_nmd_predictions:
    """
    Analyze NMD predictions from SQANTI3 classification.

    Correlates NMD status with introner presence using the original
    per-strain analysis structure: each strain's SQANTI3 file is analyzed
    against its own GTF and strain-specific introner BED, with mating-type
    region genes filtered before multi-exon filtering.

    Uses Fisher's exact test for statistical significance.
    """
    input:
        sqanti_ccmp1545 = PROJECT_ROOT / "data" / "sqanti3_output_834" / "834_isoforms_classification.filtered.txt",
        sqanti_rcc1614 = PROJECT_ROOT / "data" / "sqanti3_output_1614" / "1614_isoforms_classification.filtered.txt",
        sqanti_rcc1749 = PROJECT_ROOT / "data" / "sqanti3_output_1749" / "1749_isoforms_classification.filtered.txt",
        gtf_ccmp1545 = ANNOTATIONS_DIR / "CCMP1545.gtf",
        gtf_rcc1614 = ANNOTATIONS_DIR / "RCC1614.gtf",
        gtf_rcc1749 = ANNOTATIONS_DIR / "RCC1749.gtf",
        bed_ccmp1545 = GENOTYPING_DIR / "introner_loci" / "CCMP1545.introner_loci.bed",
        bed_rcc1614 = GENOTYPING_DIR / "introner_loci" / "RCC1614.introner_loci.bed",
        bed_rcc1749 = GENOTYPING_DIR / "introner_loci" / "RCC1749.introner_loci.bed",
        script = PROJECT_ROOT / "scripts" / "expression" / "analyze_nmd_predictions.py"
    output:
        analysis = NMD_DIR / "nmd_introner_analysis.txt",
        contingency = NMD_DIR / "nmd_contingency_tables.csv",
        plot = FIGURES_DIR / "nmd_introner_association.pdf",
        png = FIGURES_DIR / "nmd_introner_association.png"
    params:
        # Include only high-confidence structural categories
        categories = "full-splice_match,novel_in_catalog,novel_not_in_catalog"
    log:
        EXPRESSION_LOG_DIR / "nmd_analysis.log"
    conda: "../envs/isoform_analysis.yaml"
    shell:
        """
        mkdir -p {NMD_DIR}

        python {input.script} \
            --mode by_strain \
            --sqanti_ccmp1545 {input.sqanti_ccmp1545} \
            --sqanti_rcc1614 {input.sqanti_rcc1614} \
            --sqanti_rcc1749 {input.sqanti_rcc1749} \
            --gtf_ccmp1545 {input.gtf_ccmp1545} \
            --gtf_rcc1614 {input.gtf_rcc1614} \
            --gtf_rcc1749 {input.gtf_rcc1749} \
            --bed_ccmp1545 {input.bed_ccmp1545} \
            --bed_rcc1614 {input.bed_rcc1614} \
            --bed_rcc1749 {input.bed_rcc1749} \
            --mt_gtf {input.gtf_ccmp1545} \
            --output {output.analysis} \
            --contingency {output.contingency} \
            --plot {output.plot} \
            --categories {params.categories} \
            2> {log}
        """


rule nmd_by_strain:
    """
    Analyze NMD predictions stratified by strain.

    Compares NMD patterns across CCMP1545, RCC1614, and RCC1749 using
    strain-specific SQANTI3 files, GTFs, and introner BEDs.
    """
    input:
        sqanti_ccmp1545 = PROJECT_ROOT / "data" / "sqanti3_output_834" / "834_isoforms_classification.filtered.txt",
        sqanti_rcc1614 = PROJECT_ROOT / "data" / "sqanti3_output_1614" / "1614_isoforms_classification.filtered.txt",
        sqanti_rcc1749 = PROJECT_ROOT / "data" / "sqanti3_output_1749" / "1749_isoforms_classification.filtered.txt",
        gtf_ccmp1545 = ANNOTATIONS_DIR / "CCMP1545.gtf",
        gtf_rcc1614 = ANNOTATIONS_DIR / "RCC1614.gtf",
        gtf_rcc1749 = ANNOTATIONS_DIR / "RCC1749.gtf",
        bed_ccmp1545 = GENOTYPING_DIR / "introner_loci" / "CCMP1545.introner_loci.bed",
        bed_rcc1614 = GENOTYPING_DIR / "introner_loci" / "RCC1614.introner_loci.bed",
        bed_rcc1749 = GENOTYPING_DIR / "introner_loci" / "RCC1749.introner_loci.bed",
        script = PROJECT_ROOT / "scripts" / "expression" / "analyze_nmd_predictions.py"
    output:
        analysis = NMD_DIR / "nmd_by_strain_analysis.txt",
        plot = FIGURES_DIR / "nmd_by_strain.pdf",
        png = FIGURES_DIR / "nmd_by_strain.png"
    log:
        EXPRESSION_LOG_DIR / "nmd_by_strain.log"
    conda: "../envs/isoform_analysis.yaml"
    shell:
        """
        python {input.script} \
            --mode by_strain \
            --sqanti_ccmp1545 {input.sqanti_ccmp1545} \
            --sqanti_rcc1614 {input.sqanti_rcc1614} \
            --sqanti_rcc1749 {input.sqanti_rcc1749} \
            --gtf_ccmp1545 {input.gtf_ccmp1545} \
            --gtf_rcc1614 {input.gtf_rcc1614} \
            --gtf_rcc1749 {input.gtf_rcc1749} \
            --bed_ccmp1545 {input.bed_ccmp1545} \
            --bed_rcc1614 {input.bed_rcc1614} \
            --bed_rcc1749 {input.bed_rcc1749} \
            --mt_gtf {input.gtf_ccmp1545} \
            --output {output.analysis} \
            --plot {output.plot} \
            2> {log}
        """


# ============================================================
# LOCUS-AWARE ISOFORM TURNOVER DATA
# ============================================================

rule prepare_isoform_turnover_data:
    """
    Prepare the common-gene, locus-wise table used by isoform turnover models.

    SQANTI assignments are audited against sample-specific GTFs. Genotype loci
    are mapped to gene spans in each sample, missing calls remain missing, and
    callable CCMP1545-relative gains and losses are counted independently.
    """
    input:
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        ccmp1545_gtf = ANNOTATIONS_DIR / "CCMP1545.gtf",
        rcc1614_gtf = ANNOTATIONS_DIR / "RCC1614.gtf",
        rcc1749_gtf = ANNOTATIONS_DIR / "RCC1749.gtf",
        ccmp1545_sqanti = PROJECT_ROOT / "data" / "sqanti3_output_834" / "834_isoforms_classification.filtered.txt",
        rcc1614_sqanti = PROJECT_ROOT / "data" / "sqanti3_output_1614" / "1614_isoforms_classification.filtered.txt",
        rcc1749_sqanti = PROJECT_ROOT / "data" / "sqanti3_output_1749" / "1749_isoforms_classification.filtered.txt",
        script = PROJECT_ROOT / "scripts" / "expression" / "isoform_analysis" / "prepare_isoform_turnover_data.py"
    output:
        audit = ISOFORM_TURNOVER_DIR / "sqanti_gtf_audit.tsv",
        crosswalk = ISOFORM_TURNOVER_DIR / "sample_gene_crosswalk.tsv",
        richness = ISOFORM_TURNOVER_DIR / "isoform_richness_by_gene_strain.tsv",
        locus_map = ISOFORM_TURNOVER_DIR / "locus_gene_map.tsv",
        locus_events = ISOFORM_TURNOVER_DIR / "locus_level_reference_comparisons.tsv",
        introner_features = ISOFORM_TURNOVER_DIR / "introner_features_by_gene_strain.tsv",
        model_data = ISOFORM_TURNOVER_DIR / "isoform_introner_model_data.tsv",
        summary = ISOFORM_TURNOVER_DIR / "data_preparation_summary.txt"
    log:
        EXPRESSION_LOG_DIR / "prepare_isoform_turnover_data.log"
    conda: "../envs/isoform_analysis.yaml"
    shell:
        """
        mkdir -p {ISOFORM_TURNOVER_DIR}
        python {input.script} \
            --project-root {PROJECT_ROOT} \
            --genotype-matrix {input.genotype_matrix} \
            --ccmp1545-gtf {input.ccmp1545_gtf} \
            --rcc1614-gtf {input.rcc1614_gtf} \
            --rcc1749-gtf {input.rcc1749_gtf} \
            --ccmp1545-sqanti {input.ccmp1545_sqanti} \
            --rcc1614-sqanti {input.rcc1614_sqanti} \
            --rcc1749-sqanti {input.rcc1749_sqanti} \
            --mt-scaffold {TURNOVER_MT_CONTIG} \
            --mt-start {MT_START} \
            --mt-end {MT_END} \
            --output-dir {ISOFORM_TURNOVER_DIR} \
            > {log} 2>&1
        """


rule plot_introner_isoform_tracks:
    """
    Render representative sample-aligned gene model and R2C2 isoform tracks.
    """
    input:
        genotype_matrix = GENOTYPING_DIR / "genotype_matrix.final.tsv",
        introner_intron_table = RESULTS / "expression" / "functional" / "isoform_candidates" / "current_genotype_introner_vs_annotated_intron_size.tsv",
        ccmp1545_gtf = ANNOTATIONS_DIR / "CCMP1545.gtf",
        rcc1614_gtf = ANNOTATIONS_DIR / "RCC1614.gtf",
        rcc1749_gtf = ANNOTATIONS_DIR / "RCC1749.gtf",
        sqanti_data = SQANTI_DIR / "parsed_sqanti3_data.tsv",
        script = PROJECT_ROOT / "scripts" / "expression" / "isoform_analysis" / "plot_introner_isoform_tracks.py"
    output:
        pdf = FIGURES_DIR / "introner_isoform_tracks" / "{gene}.isoform_model_tracks.pdf",
        png = FIGURES_DIR / "introner_isoform_tracks" / "{gene}.isoform_model_tracks.png"
    log:
        EXPRESSION_LOG_DIR / "introner_isoform_tracks" / "{gene}.log"
    conda: "../envs/isoform_analysis.yaml"
    shell:
        """
        mkdir -p $(dirname {output.pdf})
        mkdir -p $(dirname {log})

        python {input.script} \
            --gene '{wildcards.gene}' \
            --genotype-matrix {input.genotype_matrix} \
            --annotation-dir {ANNOTATIONS_DIR} \
            --introner-intron-table {input.introner_intron_table} \
            --outdir $(dirname {output.pdf}) \
            --output-prefix '{wildcards.gene}.isoform_model_tracks' \
            --format both \
            > {log} 2>&1
        """


rule isoform_track_figures:
    """
    Target: Representative introner gene/transcript-structure track figures.
    """
    input:
        ISOFORM_TRACK_PDFS,
        ISOFORM_TRACK_PNGS


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
        ISOFORM_TURNOVER_DIR / "isoform_introner_model_data.tsv",
        ISOFORM_TURNOVER_DIR / "data_preparation_summary.txt",
        FIGURES_DIR / "shannon_diversity_distribution.pdf",
        FIGURES_DIR / "nmd_introner_association.pdf",
        ISOFORM_TRACK_PDFS,
        ISOFORM_TRACK_PNGS


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
